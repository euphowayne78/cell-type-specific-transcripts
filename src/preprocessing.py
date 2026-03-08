"""Preprocessing: normalization, QC, and basic filtering of pseudobulk data."""

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

def log1p_transform(mean_cpm: pd.DataFrame) -> pd.DataFrame:
    """Apply log1p transformation to CPM values.

    ``log1p(x) = log2(1 + x)``

    Parameters
    ----------
    mean_cpm : pd.DataFrame
        Mean CPM matrix (genes × cell types).

    Returns
    -------
    pd.DataFrame
        Log-transformed matrix (same shape / index).
    """
    result = np.log2(1 + mean_cpm)
    logger.info("Applied log1p (base-2) transform")
    return result


# ---------------------------------------------------------------------------
# Gene-level QC
# ---------------------------------------------------------------------------

def filter_low_expression(
    mean_cpm: pd.DataFrame,
    det_df: pd.DataFrame,
    min_cpm: float = 1.0,
    min_det_rate: float = 0.01,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove genes with negligible expression across all cell types.

    A gene is *kept* if in **at least one** cell type it satisfies:
    - mean CPM ≥ *min_cpm*, **and**
    - detection rate ≥ *min_det_rate*.

    Parameters
    ----------
    mean_cpm : pd.DataFrame
        Mean CPM matrix.
    det_df : pd.DataFrame
        Detection-rate matrix (same shape).
    min_cpm : float
        Minimum mean CPM in at least one cell type.
    min_det_rate : float
        Minimum detection rate in at least one cell type.

    Returns
    -------
    mean_filtered, det_filtered : pd.DataFrame
    """
    expr_pass = (mean_cpm >= min_cpm).any(axis=1)
    det_pass = (det_df >= min_det_rate).any(axis=1)
    keep = expr_pass & det_pass

    n_removed = (~keep).sum()
    logger.info(
        "Low-expression filter: kept %d / %d genes (removed %d)",
        keep.sum(), len(keep), n_removed,
    )
    return mean_cpm.loc[keep], det_df.loc[keep]


def filter_by_biotype(
    mean_df: pd.DataFrame,
    var_df: pd.DataFrame,
    allowed_biotypes: Optional[set[str]] = None,
) -> pd.DataFrame:
    """Keep only genes matching allowed biotypes.

    Parameters
    ----------
    mean_df : pd.DataFrame
        Expression matrix indexed by gene name.
    var_df : pd.DataFrame
        Gene metadata with a ``feature_biotype`` column (if available from
        Census or Ensembl BioMart annotation).
    allowed_biotypes : set[str] | None
        E.g. ``{"protein_coding", "lncRNA"}``.  *None* skips filtering.

    Returns
    -------
    pd.DataFrame
        Subset of *mean_df* matching biotypes.
    """
    if allowed_biotypes is None:
        return mean_df

    if "feature_biotype" not in var_df.columns:
        logger.warning(
            "var_df has no 'feature_biotype' column — skipping biotype filter. "
            "Annotate genes via Ensembl BioMart if biotype filtering is needed."
        )
        return mean_df

    # Map gene names to biotype
    name_to_biotype = var_df.set_index("feature_name")["feature_biotype"]
    gene_biotypes = mean_df.index.map(name_to_biotype)
    keep = gene_biotypes.isin(allowed_biotypes)

    logger.info(
        "Biotype filter (%s): kept %d / %d genes",
        ", ".join(sorted(allowed_biotypes)), keep.sum(), len(keep),
    )
    return mean_df.loc[keep.values]


# ---------------------------------------------------------------------------
# Housekeeping gene removal
# ---------------------------------------------------------------------------

# Compact housekeeping gene list (Eisenberg & Levanon 2013 + common additions).
# This is a representative subset; the full list can be loaded from a reference file.
HOUSEKEEPING_GENES: set[str] = {
    "ACTB", "GAPDH", "TUBB", "TUBA1A", "RPL13A", "RPL19", "RPS18",
    "RPS27A", "UBC", "UBB", "B2M", "HPRT1", "HMBS", "TBP", "YWHAZ",
    "SDHA", "PPIA", "RPLP0", "PGK1", "GUSB", "HSP90AB1", "LDHA",
    "NONO", "EEF1A1", "EEF2", "EIF4A2", "ATP5F1B", "HNRNPA1",
}

# Ribosomal protein genes (prefix-based)
RIBOSOMAL_PREFIXES = ("RPL", "RPS", "MRPL", "MRPS")

# Mitochondrially-encoded genes
MT_PREFIX = "MT-"


def flag_housekeeping_and_ribosomal(gene_names: pd.Index) -> pd.Series:
    """Return a boolean Series: True for housekeeping / ribosomal / MT genes.

    Parameters
    ----------
    gene_names : pd.Index
        Gene symbols.

    Returns
    -------
    pd.Series
        Boolean mask (True = housekeeping / ribosomal / mitochondrial).
    """
    is_hk = gene_names.isin(HOUSEKEEPING_GENES)
    is_ribo = gene_names.str.startswith(RIBOSOMAL_PREFIXES)
    is_mt = gene_names.str.startswith(MT_PREFIX)
    return is_hk | is_ribo | is_mt


def remove_housekeeping(
    mean_df: pd.DataFrame,
    det_df: pd.DataFrame,
    remove_ribosomal: bool = True,
    remove_mitochondrial: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove housekeeping, ribosomal, and mitochondrial genes.

    Parameters
    ----------
    mean_df, det_df : pd.DataFrame
        Expression and detection-rate matrices.
    remove_ribosomal : bool
        Also remove ribosomal protein genes (RPL*, RPS*, MRPL*, MRPS*).
    remove_mitochondrial : bool
        Also remove MT-* genes.

    Returns
    -------
    mean_filtered, det_filtered : pd.DataFrame
    """
    mask = mean_df.index.isin(HOUSEKEEPING_GENES)
    if remove_ribosomal:
        mask = mask | mean_df.index.str.startswith(RIBOSOMAL_PREFIXES)
    if remove_mitochondrial:
        mask = mask | mean_df.index.str.startswith(MT_PREFIX)

    n_removed = mask.sum()
    logger.info(
        "Housekeeping filter: removed %d genes (HK=%d, ribo=%s, MT=%s)",
        n_removed,
        mean_df.index.isin(HOUSEKEEPING_GENES).sum(),
        remove_ribosomal,
        remove_mitochondrial,
    )
    return mean_df.loc[~mask], det_df.loc[~mask]


# ---------------------------------------------------------------------------
# Sex-chromosome filter
# ---------------------------------------------------------------------------

# Genes on chrX / chrY that may confound specificity analysis.
# A proper implementation would load a chromosome map from Ensembl BioMart.
# For now, provide a utility that accepts an external gene→chrom mapping.

def remove_sex_chromosome_genes(
    mean_df: pd.DataFrame,
    det_df: pd.DataFrame,
    gene_chrom_map: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Remove genes on chrX and chrY.

    Parameters
    ----------
    gene_chrom_map : pd.Series
        Series mapping gene name → chromosome (e.g. ``"1"``, ``"X"``).
    """
    chroms = mean_df.index.map(gene_chrom_map)
    keep = ~chroms.isin({"X", "Y"})
    n_removed = (~keep).sum()
    logger.info("Sex-chromosome filter: removed %d genes", n_removed)
    return mean_df.loc[keep.values], det_df.loc[keep.values]


# ---------------------------------------------------------------------------
# Pipeline convenience
# ---------------------------------------------------------------------------

def run_preprocessing(
    mean_df: pd.DataFrame,
    det_df: pd.DataFrame,
    min_cpm: float = 1.0,
    min_det_rate: float = 0.01,
    remove_hk: bool = True,
    log_transform: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the standard preprocessing pipeline.

    Steps:
    1. Filter low-expression genes.
    2. Remove housekeeping / ribosomal / MT genes.
    3. (Optional) log1p-transform mean expression.

    Returns
    -------
    mean_processed, det_processed : pd.DataFrame
    """
    mean_out, det_out = filter_low_expression(
        mean_df, det_df, min_cpm=min_cpm, min_det_rate=min_det_rate,
    )
    if remove_hk:
        mean_out, det_out = remove_housekeeping(mean_out, det_out)
    if log_transform:
        mean_out = log1p_transform(mean_out)

    logger.info(
        "Preprocessing complete: %d genes × %d cell types",
        mean_out.shape[0], mean_out.shape[1],
    )
    return mean_out, det_out
