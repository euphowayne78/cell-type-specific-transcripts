"""CELLxGENE Census data acquisition and pseudobulk aggregation."""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default query filters
# ---------------------------------------------------------------------------
DEFAULT_OBS_VALUE_FILTER = (
    "disease == 'normal' "
    "and is_primary_data == True"
)

# Patterns in development_stage to exclude (fetal/embryonic/pediatric)
EXCLUDE_DEV_STAGE_PATTERNS = (
    "embryonic|fetal|Carnegie|gestational|prenatal|"
    "newborn|infant|child|adolescent"
)

EXCLUDE_TISSUE_GENERAL = {"embryo", "organoid", "cell culture"}

OBS_COLUMNS = [
    "soma_joinid",
    "cell_type",
    "cell_type_ontology_term_id",
    "tissue",
    "tissue_general",
    "disease",
    "development_stage",
    "sex",
    "dataset_id",
    "assay",
    "is_primary_data",
]

MIN_CELLS_PER_TYPE = 500
MIN_DATASETS_PER_TYPE = 2
MAX_CELLS_PER_TYPE = 10_000  # subsample cap for memory


# ---------------------------------------------------------------------------
# Census helpers
# ---------------------------------------------------------------------------

def get_census_version() -> str:
    """Return the latest stable Census version string."""
    import cellxgene_census
    desc = cellxgene_census.get_census_version_description("stable")
    # API returns different keys depending on version - try common ones
    for key in ["census_version", "release_build", "alias"]:
        if key in desc:
            return str(desc[key])
    # Fallback: return the whole dict as string
    return str(desc)


def open_census(census_version: str = "stable"):
    """Open and return a Census SOMA collection.

    Use as a context manager::

        with open_census() as census:
            ...
    """
    import cellxgene_census
    return cellxgene_census.open_soma(census_version=census_version)


# ---------------------------------------------------------------------------
# Observation metadata
# ---------------------------------------------------------------------------

def get_obs_metadata(
    census,
    obs_value_filter: str = DEFAULT_OBS_VALUE_FILTER,
    columns: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Query cell-level metadata from Census.

    Parameters
    ----------
    census : soma.Collection
        An open Census handle (from ``open_census``).
    obs_value_filter : str
        SOMA value-filter expression applied server-side.
    columns : list[str] | None
        Columns to retrieve.  ``None`` ⇒ :data:`OBS_COLUMNS`.

    Returns
    -------
    pd.DataFrame
        One row per cell passing the filter.
    """
    if columns is None:
        columns = OBS_COLUMNS

    logger.info("Querying obs metadata with filter: %s", obs_value_filter)
    obs_df = (
        census["census_data"]["homo_sapiens"]
        .obs.read(value_filter=obs_value_filter, column_names=columns)
        .concat()
        .to_pandas()
    )
    logger.info("Retrieved %s cells", f"{len(obs_df):,}")
    return obs_df


def get_var_metadata(census) -> pd.DataFrame:
    """Retrieve gene (var) metadata from Census.

    Returns
    -------
    pd.DataFrame
        Indexed by ``soma_joinid`` with columns ``feature_id``,
        ``feature_name``, ``feature_length``.
    """
    var_df = (
        census["census_data"]["homo_sapiens"]
        .ms["RNA"]
        .var.read(
            column_names=["soma_joinid", "feature_id", "feature_name", "feature_length"]
        )
        .concat()
        .to_pandas()
    )
    var_df = var_df.set_index("soma_joinid")
    logger.info("Retrieved %s genes", f"{len(var_df):,}")
    return var_df


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------

def filter_adult_cells(obs_df: pd.DataFrame) -> pd.DataFrame:
    """Remove fetal / embryonic / pediatric cells and excluded tissues.

    Parameters
    ----------
    obs_df : pd.DataFrame
        Output of :func:`get_obs_metadata`.

    Returns
    -------
    pd.DataFrame
        Filtered copy.
    """
    dev_mask = obs_df["development_stage"].str.contains(
        EXCLUDE_DEV_STAGE_PATTERNS, case=False, na=False,
    )
    tissue_mask = obs_df["tissue_general"].str.lower().isin(
        {t.lower() for t in EXCLUDE_TISSUE_GENERAL}
    )
    adult_df = obs_df[~dev_mask & ~tissue_mask].copy()
    n_removed = len(obs_df) - len(adult_df)
    logger.info(
        "Filtered to %s adult cells (removed %s fetal/excluded)",
        f"{len(adult_df):,}", f"{n_removed:,}",
    )
    return adult_df


def get_cell_type_summary(
    obs_df: pd.DataFrame,
    min_cells: int = MIN_CELLS_PER_TYPE,
    min_datasets: int = MIN_DATASETS_PER_TYPE,
) -> pd.DataFrame:
    """Summarise cell-type representation and apply quality thresholds.

    Parameters
    ----------
    obs_df : pd.DataFrame
        Filtered observation metadata.
    min_cells : int
        Minimum number of cells required per cell type.
    min_datasets : int
        Minimum number of independent datasets required.

    Returns
    -------
    pd.DataFrame
        Columns: ``cell_type``, ``cell_type_ontology_term_id``,
        ``n_cells``, ``n_datasets``, ``tissues``.  Sorted by *n_cells*
        descending.
    """
    summary = (
        obs_df
        .groupby(["cell_type", "cell_type_ontology_term_id"])
        .agg(
            n_cells=("soma_joinid", "count"),
            n_datasets=("dataset_id", "nunique"),
            tissues=("tissue_general", lambda s: sorted(s.unique())),
        )
        .reset_index()
    )
    filtered = summary[
        (summary["n_cells"] >= min_cells)
        & (summary["n_datasets"] >= min_datasets)
    ].sort_values("n_cells", ascending=False).reset_index(drop=True)

    logger.info(
        "Cell types passing thresholds (>=%d cells, >=%d datasets): %d / %d",
        min_cells, min_datasets, len(filtered), len(summary),
    )
    return filtered


# ---------------------------------------------------------------------------
# Pseudobulk computation
# ---------------------------------------------------------------------------

def compute_pseudobulk(
    census,
    obs_df: pd.DataFrame,
    cell_types: list[str],
    max_cells_per_type: int = MAX_CELLS_PER_TYPE,
    random_seed: int = 42,
    checkpoint_dir: Optional[Path] = None,
    checkpoint_every: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute per-cell-type pseudobulk expression and detection rate.

    For each cell type the function:

    1. Retrieves the raw count matrix (subsampled if necessary).
    2. CPM-normalises each cell independently.
    3. Computes the **mean CPM** across cells  → *mean expression*.
    4. Computes the **fraction of cells with count > 0** → *detection rate*.

    Parameters
    ----------
    census : soma.Collection
        Open Census handle.
    obs_df : pd.DataFrame
        Filtered observation metadata (only adult / healthy cells).
    cell_types : list[str]
        Ordered list of cell-type labels to process.
    max_cells_per_type : int
        If a cell type has more cells than this value, randomly subsample.
    random_seed : int
        Seed for reproducible subsampling.

    Returns
    -------
    mean_df : pd.DataFrame
        Shape ``(n_genes, n_cell_types)``, index = gene symbol.
        Values are mean CPM.
    det_df : pd.DataFrame
        Same shape.  Values in [0, 1] — fraction of cells expressing.
    """
    import cellxgene_census
    from scipy.sparse import issparse

    rng = np.random.default_rng(random_seed)

    mean_expressions: dict[str, np.ndarray] = {}
    detection_rates: dict[str, np.ndarray] = {}
    gene_names: Optional[np.ndarray] = None
    gene_ids: Optional[np.ndarray] = None
    start_idx = 1

    # Load from checkpoint if available
    if checkpoint_dir:
        checkpoint_dir = Path(checkpoint_dir)
        checkpoint_files = sorted(checkpoint_dir.glob("checkpoint_*.pkl"))
        if checkpoint_files:
            latest_checkpoint = checkpoint_files[-1]
            logger.info("Resuming from checkpoint: %s", latest_checkpoint)
            import pickle
            with open(latest_checkpoint, "rb") as f:
                ckpt = pickle.load(f)
                mean_expressions = ckpt["mean_expressions"]
                detection_rates = ckpt["detection_rates"]
                gene_names = ckpt["gene_names"]
                gene_ids = ckpt["gene_ids"]
                start_idx = ckpt["processed_idx"] + 1
            logger.info("Resuming from cell type %d/%d", start_idx, len(cell_types))

    for idx, cell_type in enumerate(cell_types, 1):
        if idx < start_idx:
            continue  # Skip already processed
        logger.info("[%d/%d] Processing: %s", idx, len(cell_types), cell_type)

        cell_ids = obs_df.loc[
            obs_df["cell_type"] == cell_type, "soma_joinid"
        ].values

        if len(cell_ids) == 0:
            logger.warning("  No cells found — skipping.")
            continue

        # Subsample if needed
        if len(cell_ids) > max_cells_per_type:
            cell_ids = rng.choice(cell_ids, size=max_cells_per_type, replace=False)
            logger.info("  Subsampled to %d cells", max_cells_per_type)

        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_coords=cell_ids.tolist(),
            column_names={
                "obs": ["soma_joinid", "cell_type"],
                "var": ["soma_joinid", "feature_id", "feature_name"],
            },
        )

        X = adata.X
        if issparse(X):
            cell_totals = np.array(X.sum(axis=1)).ravel()
            cell_totals[cell_totals == 0] = 1
            # Sparse CPM
            X_cpm = X.multiply(1e6 / cell_totals[:, np.newaxis])
            mean_expr = np.asarray(X_cpm.mean(axis=0)).ravel()
            det_rate = np.asarray((X > 0).mean(axis=0)).ravel()
        else:
            cell_totals = X.sum(axis=1, keepdims=True)
            cell_totals[cell_totals == 0] = 1
            X_cpm = X * (1e6 / cell_totals)
            mean_expr = X_cpm.mean(axis=0)
            det_rate = (X > 0).astype(float).mean(axis=0)

        mean_expressions[cell_type] = mean_expr
        detection_rates[cell_type] = det_rate

        if gene_names is None:
            gene_names = adata.var["feature_name"].values
            gene_ids = adata.var["feature_id"].values

        logger.info(
            "  Done — %d cells, %d genes detected (>0 in ≥1 cell)",
            len(cell_ids), int((det_rate > 0).sum()),
        )

        # Save checkpoint periodically
        if checkpoint_dir and idx % checkpoint_every == 0:
            checkpoint_path = Path(checkpoint_dir) / f"checkpoint_{idx}.pkl"
            checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
            import pickle
            with open(checkpoint_path, "wb") as f:
                pickle.dump({
                    "mean_expressions": mean_expressions,
                    "detection_rates": detection_rates,
                    "gene_names": gene_names,
                    "gene_ids": gene_ids,
                    "processed_idx": idx,
                    "cell_types": cell_types[:idx],
                }, f)
            logger.info("  ✓ Checkpoint saved: %s", checkpoint_path)

    # Assemble DataFrames
    mean_df = pd.DataFrame(mean_expressions, index=gene_names)
    det_df = pd.DataFrame(detection_rates, index=gene_names)
    mean_df.index.name = "gene_name"
    det_df.index.name = "gene_name"

    # Store Ensembl IDs as an attribute (accessible via df.attrs)
    id_map = pd.Series(gene_ids, index=gene_names, name="feature_id")
    mean_df.attrs["gene_id_map"] = id_map
    det_df.attrs["gene_id_map"] = id_map

    return mean_df, det_df


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def save_pseudobulk(
    mean_df: pd.DataFrame,
    det_df: pd.DataFrame,
    output_dir: Path | str,
    obs_df: Optional[pd.DataFrame] = None,
    cell_type_summary: Optional[pd.DataFrame] = None,
) -> None:
    """Save pseudobulk matrices and optional metadata to *output_dir*.

    Files written:
    - ``pseudobulk_mean_cpm.parquet``
    - ``pseudobulk_detection_rate.parquet``
    - ``obs_metadata.parquet`` (if *obs_df* provided)
    - ``cell_type_summary.parquet`` (if provided)
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Save gene_id_map separately — parquet can't serialize a Series in attrs
    gene_id_map = mean_df.attrs.get("gene_id_map")
    if gene_id_map is not None:
        gene_id_map.to_frame().to_parquet(out / "gene_id_map.parquet")

    mean_df_save = mean_df.copy()
    mean_df_save.attrs = {}
    det_df_save = det_df.copy()
    det_df_save.attrs = {}

    mean_df_save.to_parquet(out / "pseudobulk_mean_cpm.parquet")
    det_df_save.to_parquet(out / "pseudobulk_detection_rate.parquet")
    logger.info("Saved pseudobulk matrices to %s", out)

    if obs_df is not None:
        obs_df.to_parquet(out / "obs_metadata.parquet", index=False)
    if cell_type_summary is not None:
        cell_type_summary.to_parquet(out / "cell_type_summary.parquet", index=False)


def load_pseudobulk(
    data_dir: Path | str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load cached pseudobulk matrices.

    Returns
    -------
    mean_df, det_df : pd.DataFrame
        Mean CPM and detection-rate matrices.
    """
    d = Path(data_dir)
    mean_df = pd.read_parquet(d / "pseudobulk_mean_cpm.parquet")
    det_df = pd.read_parquet(d / "pseudobulk_detection_rate.parquet")
    logger.info(
        "Loaded pseudobulk: %d genes × %d cell types",
        mean_df.shape[0], mean_df.shape[1],
    )
    return mean_df, det_df
