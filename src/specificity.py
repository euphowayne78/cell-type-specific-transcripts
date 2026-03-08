"""Specificity metrics for identifying cell-type-specific transcripts."""

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Tau specificity index (Yanai et al., 2005)
# ---------------------------------------------------------------------------

def compute_tau(expression: pd.DataFrame) -> pd.Series:
    """Compute the Tau tissue/cell-type specificity index.

    .. math::

        \\tau = \\frac{\\sum_{i=1}^{n}(1 - \\hat{x}_i)}{n - 1}

    where :math:`\\hat{x}_i = x_i / \\max(x)`.

    Parameters
    ----------
    expression : pd.DataFrame
        Expression matrix (genes × cell types).  Values should be
        non-negative (e.g. mean CPM or log-CPM).

    Returns
    -------
    pd.Series
        Tau value per gene.  Range [0, 1]:
        0 = ubiquitous, 1 = perfectly cell-type specific.
    """
    max_expr = expression.max(axis=1)
    n = expression.shape[1]

    # Genes with zero expression everywhere → Tau = 0
    zero_mask = max_expr == 0
    max_safe = max_expr.copy()
    max_safe[zero_mask] = 1  # prevent division by zero

    x_hat = expression.div(max_safe, axis=0)
    tau = (1.0 - x_hat).sum(axis=1) / (n - 1)
    tau[zero_mask] = 0.0  # explicitly set zero-expression genes to Tau = 0
    tau.name = "tau"
    return tau


# ---------------------------------------------------------------------------
# Fold-change vs. next-highest cell type
# ---------------------------------------------------------------------------

def compute_fold_change(
    expression: pd.DataFrame,
    pseudocount: float = 0.1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute fold-change of each gene in each cell type vs. the next highest.

    Parameters
    ----------
    expression : pd.DataFrame
        Expression matrix (genes × cell types).
    pseudocount : float
        Added to both numerator and denominator to avoid division by zero
        and to dampen fold-change for very low expression values.

    Returns
    -------
    fc_df : pd.DataFrame
        Same shape as *expression*.  Value = fold-change over next-highest.
    second_df : pd.DataFrame
        Expression of the second-highest cell type for each gene.
    """
    fc_data = {}
    second_data = {}

    for ct in expression.columns:
        target = expression[ct]
        others_max = expression.drop(columns=ct).max(axis=1)
        fc = (target + pseudocount) / (others_max + pseudocount)
        fc_data[ct] = fc
        second_data[ct] = others_max

    fc_df = pd.DataFrame(fc_data, index=expression.index)
    second_df = pd.DataFrame(second_data, index=expression.index)
    return fc_df, second_df


# ---------------------------------------------------------------------------
# Detection rate thresholds
# ---------------------------------------------------------------------------

def compute_detection_specificity(
    det_df: pd.DataFrame,
    on_target_min: float = 0.25,
    off_target_max: float = 0.05,
) -> pd.DataFrame:
    """Flag genes passing detection-rate thresholds per cell type.

    For each (gene, cell_type) pair, checks:
    - detection rate in *cell_type* ≥ ``on_target_min``
    - detection rate in **all other** cell types ≤ ``off_target_max``

    Parameters
    ----------
    det_df : pd.DataFrame
        Detection-rate matrix (genes × cell types), values in [0, 1].

    Returns
    -------
    pd.DataFrame
        Boolean matrix — True where both thresholds are met.
    """
    result = pd.DataFrame(False, index=det_df.index, columns=det_df.columns)

    for ct in det_df.columns:
        on_target = det_df[ct] >= on_target_min
        off_target = det_df.drop(columns=ct).max(axis=1) <= off_target_max
        result[ct] = on_target & off_target

    return result


# ---------------------------------------------------------------------------
# AUROC (one-vs-rest)
# ---------------------------------------------------------------------------

def compute_auroc_per_gene(
    expression: pd.DataFrame,
    det_df: pd.DataFrame,
) -> pd.DataFrame:
    """Compute one-vs-rest AUROC for each gene × cell type.

    Uses the pseudobulk expression as a classifier score.  AUROC measures
    how well a single gene's expression separates one cell type from all
    others.

    Since we operate on pseudobulk (one value per cell type, not per cell),
    this is an approximation.  For rigorous per-cell AUROC, use the raw
    single-cell matrix.

    Parameters
    ----------
    expression : pd.DataFrame
        Expression matrix (genes × cell types).
    det_df : pd.DataFrame
        Detection-rate matrix (used as secondary signal).

    Returns
    -------
    pd.DataFrame
        AUROC values, shape (n_genes, n_cell_types).
        Values in [0, 1]; 0.5 = random, 1.0 = perfect separation.
    """
    from sklearn.metrics import roc_auc_score

    n_types = expression.shape[1]
    if n_types < 3:
        logger.warning("AUROC requires ≥3 cell types; returning NaN")
        return pd.DataFrame(
            np.nan, index=expression.index, columns=expression.columns,
        )

    auroc_data = {}
    for ct in expression.columns:
        y_true = np.zeros(n_types)
        y_true[list(expression.columns).index(ct)] = 1

        gene_aurocs = []
        for gene_idx in range(len(expression)):
            y_score = expression.iloc[gene_idx].values
            if y_score.max() == y_score.min():
                gene_aurocs.append(0.5)
            else:
                try:
                    auc = roc_auc_score(y_true, y_score)
                    gene_aurocs.append(auc)
                except ValueError:
                    gene_aurocs.append(0.5)
        auroc_data[ct] = gene_aurocs

    return pd.DataFrame(auroc_data, index=expression.index)


# ---------------------------------------------------------------------------
# Combined specificity scoring
# ---------------------------------------------------------------------------

def score_specificity(
    mean_cpm: pd.DataFrame,
    det_df: pd.DataFrame,
    tau_threshold: float = 0.85,
    fc_threshold: float = 10.0,
    on_target_det_min: float = 0.25,
    off_target_det_max: float = 0.05,
    on_target_expr_min: float = 1.0,
    pseudocount: float = 0.1,
) -> pd.DataFrame:
    """Run all specificity metrics and return a long-form results table.

    Each row is a (gene, cell_type) pair that passes the Tau threshold.
    Additional columns flag whether each supporting metric is met.

    Parameters
    ----------
    mean_cpm : pd.DataFrame
        Mean CPM expression (genes × cell types).  **Not** log-transformed
        for fold-change computation; log values used for Tau.
    det_df : pd.DataFrame
        Detection-rate matrix.

    Returns
    -------
    pd.DataFrame
        Columns: ``gene``, ``cell_type``, ``tau``, ``mean_cpm``,
        ``detection_rate``, ``fold_change``, ``second_highest_cpm``,
        ``on_target_expr_pass``, ``on_target_det_pass``,
        ``off_target_det_pass``, ``fc_pass``, ``all_pass``.
    """
    # Deduplicate gene index — keep max expression row per duplicate gene name
    if mean_cpm.index.duplicated().any():
        n_dups = mean_cpm.index.duplicated().sum()
        logger.warning("Duplicate gene names detected (%d); keeping max-expression row per gene.", n_dups)
        mean_cpm = mean_cpm.groupby(level=0).max()
        det_df = det_df.groupby(level=0).max()

    # Log-transform for Tau (Tau works best on log-scale)
    log_expr = np.log2(1 + mean_cpm)

    # --- Tau ---
    tau = compute_tau(log_expr)

    # --- Fold-change (on raw CPM for biological interpretability) ---
    fc_df, second_df = compute_fold_change(mean_cpm, pseudocount=pseudocount)

    # --- Detection specificity ---
    det_pass_matrix = compute_detection_specificity(
        det_df, on_target_min=on_target_det_min, off_target_max=off_target_det_max,
    )

    # --- Collect results in long form ---
    # For each gene, identify the cell type with max expression
    best_ct = mean_cpm.idxmax(axis=1)

    records = []
    for gene in mean_cpm.index:
        t = tau.loc[gene]
        if t < tau_threshold:
            continue  # Skip non-specific genes

        ct = best_ct.loc[gene]
        expr_val = mean_cpm.loc[gene, ct]
        det_val = det_df.loc[gene, ct]
        fc_val = fc_df.loc[gene, ct]
        second_val = second_df.loc[gene, ct]
        max_off_det = det_df.loc[gene].drop(ct).max()

        records.append({
            "gene": gene,
            "cell_type": ct,
            "tau": round(t, 4),
            "mean_cpm": round(expr_val, 2),
            "log2_cpm": round(np.log2(1 + expr_val), 2),
            "detection_rate": round(det_val, 4),
            "fold_change": round(fc_val, 2),
            "second_highest_cpm": round(second_val, 2),
            "max_off_target_det": round(max_off_det, 4),
            "on_target_expr_pass": expr_val >= on_target_expr_min,
            "on_target_det_pass": det_val >= on_target_det_min,
            "off_target_det_pass": max_off_det <= off_target_det_max,
            "fc_pass": fc_val >= fc_threshold,
        })

    if not records:
        logger.warning("No genes passed Tau threshold %.2f", tau_threshold)
        return pd.DataFrame()

    results = pd.DataFrame(records)
    results["all_pass"] = (
        results["on_target_expr_pass"]
        & results["on_target_det_pass"]
        & results["off_target_det_pass"]
        & results["fc_pass"]
    )
    results = results.sort_values(
        ["all_pass", "tau", "fold_change"], ascending=[False, False, False],
    ).reset_index(drop=True)

    n_pass = results["all_pass"].sum()
    logger.info(
        "Specificity scoring: %d genes passed Tau > %.2f; %d passed ALL thresholds",
        len(results), tau_threshold, n_pass,
    )
    return results


# ---------------------------------------------------------------------------
# Composite specificity score
# ---------------------------------------------------------------------------

def compute_composite_score(
    results: pd.DataFrame,
    weight_tau: float = 0.35,
    weight_fc: float = 0.25,
    weight_detection: float = 0.20,
    weight_expression: float = 0.10,
    weight_off_target: float = 0.10,
) -> pd.DataFrame:
    """Add a weighted composite specificity score to the results DataFrame.

    Each metric is normalised to [0, 1] and combined as a weighted sum,
    giving a single score that can be used to rank candidates within and
    across cell types.

    Normalisation scheme:
    - ``tau``               — already in [0, 1]; used directly.
    - ``fold_change``       — log10-scaled then divided by the dataset max.
    - ``detection_rate``    — already in [0, 1]; used directly.
    - ``log2_cpm``          — divided by the dataset max.
    - ``max_off_target_det``— inverted (1 − value) so lower off-target → higher score.

    Parameters
    ----------
    results : pd.DataFrame
        Output of :func:`score_specificity`.
    weight_tau : float
        Weight for Tau (default 0.35).
    weight_fc : float
        Weight for fold-change (default 0.25).
    weight_detection : float
        Weight for on-target detection rate (default 0.20).
    weight_expression : float
        Weight for on-target log2 CPM (default 0.10).
    weight_off_target : float
        Weight for inverted off-target detection rate (default 0.10).

    Returns
    -------
    pd.DataFrame
        Input DataFrame with an additional ``composite_score`` column [0, 1],
        re-sorted by ``all_pass`` (True first) then ``composite_score`` descending.
    """
    if results.empty:
        return results

    df = results.copy()

    # Tau: already [0, 1]
    tau_norm = df["tau"]

    # Fold-change: log10 scale to handle wide dynamic range, then normalise
    fc_log = np.log10(df["fold_change"].clip(lower=1.0) + 1.0)
    fc_max = fc_log.max()
    fc_norm = fc_log / fc_max if fc_max > 0 else fc_log

    # On-target detection rate: already [0, 1]
    det_norm = df["detection_rate"]

    # On-target expression: normalise log2_cpm to [0, 1]
    expr_max = df["log2_cpm"].max()
    expr_norm = df["log2_cpm"] / expr_max if expr_max > 0 else df["log2_cpm"]

    # Off-target detection: invert so lower off-target → higher contribution
    off_norm = 1.0 - df["max_off_target_det"]

    df["composite_score"] = (
        weight_tau * tau_norm
        + weight_fc * fc_norm
        + weight_detection * det_norm
        + weight_expression * expr_norm
        + weight_off_target * off_norm
    ).round(4)

    df = df.sort_values(
        ["all_pass", "composite_score"], ascending=[False, False]
    ).reset_index(drop=True)

    logger.info(
        "Composite scores computed. Range: %.3f – %.3f",
        df["composite_score"].min(), df["composite_score"].max(),
    )
    return df


# ---------------------------------------------------------------------------
# Per-cell-type query helper
# ---------------------------------------------------------------------------

def query_cell_type(
    results: pd.DataFrame,
    cell_type: str,
    all_pass_only: bool = True,
    top_n: int = 50,
) -> pd.DataFrame:
    """Return top specific transcripts for a given cell type.

    Parameters
    ----------
    results : pd.DataFrame
        Output of :func:`score_specificity`.
    cell_type : str
        Cell type to query.
    all_pass_only : bool
        If True, only return genes passing all thresholds.
    top_n : int
        Maximum number of genes to return.

    Returns
    -------
    pd.DataFrame
        Filtered and ranked subset.
    """
    subset = results[results["cell_type"] == cell_type].copy()
    if all_pass_only:
        subset = subset[subset["all_pass"]]
    sort_col = "composite_score" if "composite_score" in subset.columns else "tau"
    return subset.sort_values(sort_col, ascending=False).head(top_n).reset_index(drop=True)
