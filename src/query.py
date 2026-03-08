"""End-user query interface logic for cell-type specific transcript lookup."""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result loading
# ---------------------------------------------------------------------------

def load_results(results_dir: Path | str) -> pd.DataFrame:
    """Load the full specificity results table.

    Parameters
    ----------
    results_dir : Path | str
        Directory containing ``specificity_scores.csv``.

    Returns
    -------
    pd.DataFrame
    """
    path = Path(results_dir) / "specificity_scores.csv"
    df = pd.read_csv(path)
    logger.info("Loaded %d results from %s", len(df), path)
    return df


def load_expression_matrices(
    data_dir: Path | str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load processed pseudobulk matrices.

    Returns
    -------
    mean_df, det_df : pd.DataFrame
    """
    d = Path(data_dir)
    mean_df = pd.read_parquet(d / "pseudobulk_mean_cpm.parquet")
    det_df = pd.read_parquet(d / "pseudobulk_detection_rate.parquet")
    return mean_df, det_df


# ---------------------------------------------------------------------------
# Query functions
# ---------------------------------------------------------------------------

def get_available_cell_types(results: pd.DataFrame) -> list[str]:
    """Return sorted list of cell types with at least one specific transcript.

    Parameters
    ----------
    results : pd.DataFrame
        Full specificity results.

    Returns
    -------
    list[str]
    """
    return sorted(results["cell_type"].unique())


def query_transcripts(
    results: pd.DataFrame,
    cell_type: str,
    tau_min: float = 0.85,
    fc_min: float = 10.0,
    det_min: float = 0.25,
    off_target_max: float = 0.05,
    expr_min: float = 1.0,
    annotation_tier: Optional[int] = None,
    all_pass_only: bool = True,
    top_n: int = 50,
) -> pd.DataFrame:
    """Query specific transcripts for a cell type with adjustable thresholds.

    Parameters
    ----------
    results : pd.DataFrame
        Full specificity results (from ``specificity.score_specificity``
        + ``filters.apply_post_filters``).
    cell_type : str
        Target cell type.
    tau_min : float
        Minimum Tau specificity index.
    fc_min : float
        Minimum fold-change vs. next-highest cell type.
    det_min : float
        Minimum on-target detection rate.
    off_target_max : float
        Maximum off-target detection rate in any other cell type.
    expr_min : float
        Minimum on-target mean CPM.
    annotation_tier : int | None
        If set, only return genes at this tier or better (lower number = better).
    all_pass_only : bool
        If True, only genes passing all thresholds.
    top_n : int
        Max results to return.

    Returns
    -------
    pd.DataFrame
        Ranked transcripts for the queried cell type.
    """
    sub = results[results["cell_type"] == cell_type].copy()

    if sub.empty:
        logger.warning("No results for cell type: %s", cell_type)
        return pd.DataFrame()

    # Apply threshold filters
    mask = (
        (sub["tau"] >= tau_min)
        & (sub["fold_change"] >= fc_min)
        & (sub["detection_rate"] >= det_min)
        & (sub["max_off_target_det"] <= off_target_max)
        & (sub["mean_cpm"] >= expr_min)
    )

    if annotation_tier is not None and "annotation_tier" in sub.columns:
        mask = mask & (sub["annotation_tier"] <= annotation_tier)

    filtered = sub[mask].sort_values(
        ["tau", "fold_change"], ascending=[False, False]
    ).head(top_n).reset_index(drop=True)

    logger.info(
        "Query '%s': %d transcripts (from %d candidates)",
        cell_type, len(filtered), len(sub),
    )
    return filtered


def compare_across_cell_types(
    results: pd.DataFrame,
    gene: str,
    mean_df: pd.DataFrame,
    det_df: pd.DataFrame,
) -> pd.DataFrame:
    """Show expression of a single gene across all cell types.

    Useful for verifying that a candidate transcript is truly specific.

    Parameters
    ----------
    gene : str
        Gene symbol.
    mean_df : pd.DataFrame
        Mean CPM matrix.
    det_df : pd.DataFrame
        Detection-rate matrix.

    Returns
    -------
    pd.DataFrame
        One row per cell type, sorted by expression descending.
        Columns: ``cell_type``, ``mean_cpm``, ``detection_rate``.
    """
    if gene not in mean_df.index:
        raise ValueError(f"Gene '{gene}' not found in expression data")

    comparison = pd.DataFrame({
        "cell_type": mean_df.columns,
        "mean_cpm": mean_df.loc[gene].values,
        "detection_rate": det_df.loc[gene].values,
    })
    comparison["log2_cpm"] = pd.Series(
        __import__("numpy").log2(1 + comparison["mean_cpm"])
    )
    return comparison.sort_values("mean_cpm", ascending=False).reset_index(drop=True)


def generate_report(
    results: pd.DataFrame,
    cell_type: str,
    mean_df: pd.DataFrame,
    det_df: pd.DataFrame,
    top_n: int = 20,
) -> str:
    """Generate a text summary report for a cell type.

    Parameters
    ----------
    results : pd.DataFrame
        Full specificity results.
    cell_type : str
        Target cell type.
    mean_df, det_df : pd.DataFrame
        Expression and detection-rate matrices.
    top_n : int
        Number of top hits to include.

    Returns
    -------
    str
        Formatted report text.
    """
    hits = query_transcripts(results, cell_type, top_n=top_n)

    lines = [
        f"{'='*60}",
        f"Cell-Type Specific Transcript Report: {cell_type}",
        f"{'='*60}",
        f"Total candidates passing all thresholds: {len(hits)}",
        "",
    ]

    if hits.empty:
        lines.append("No transcripts passed all specificity thresholds.")
        return "\n".join(lines)

    # Tier summary
    if "annotation_tier" in hits.columns:
        tier_counts = hits["annotation_tier"].value_counts().sort_index()
        lines.append("Annotation Tier Breakdown:")
        for tier, count in tier_counts.items():
            lines.append(f"  Tier {tier}: {count} genes")
        lines.append("")

    lines.append(f"Top {min(top_n, len(hits))} Hits:")
    lines.append(f"{'Gene':<15} {'Tau':>6} {'CPM':>10} {'Det%':>6} {'FC':>8} {'Tier':>5}")
    lines.append("-" * 55)

    for _, row in hits.head(top_n).iterrows():
        tier = row.get("annotation_tier", "?")
        lines.append(
            f"{row['gene']:<15} {row['tau']:>6.3f} {row['mean_cpm']:>10.1f} "
            f"{row['detection_rate']:>5.1%} {row['fold_change']:>8.1f} {tier:>5}"
        )

    return "\n".join(lines)
