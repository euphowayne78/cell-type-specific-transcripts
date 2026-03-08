"""Visualization: heatmaps, dot plots, and expression profiles."""

import logging
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

# Consistent style
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
})


# ---------------------------------------------------------------------------
# Expression heatmap
# ---------------------------------------------------------------------------

def plot_expression_heatmap(
    expression: pd.DataFrame,
    genes: Optional[list[str]] = None,
    cell_types: Optional[list[str]] = None,
    title: str = "Cell-Type Specific Transcript Expression",
    figsize: Optional[tuple[float, float]] = None,
    cmap: str = "YlOrRd",
    z_score: Optional[int] = 0,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot a clustered heatmap of expression (genes × cell types).

    Parameters
    ----------
    expression : pd.DataFrame
        Expression matrix (genes × cell types).
    genes : list[str] | None
        Subset of genes to display.  *None* → all.
    cell_types : list[str] | None
        Subset of cell types.  *None* → all.
    title : str
        Plot title.
    figsize : tuple | None
        Figure size.  Auto-computed if *None*.
    cmap : str
        Colormap name.
    z_score : int | None
        Apply z-score normalisation across axis (0 = per gene, 1 = per
        cell type, None = no z-score).
    save_path : str | None
        If provided, save figure to this path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    data = expression.copy()
    if genes is not None:
        data = data.loc[data.index.isin(genes)]
    if cell_types is not None:
        data = data[[c for c in cell_types if c in data.columns]]

    if figsize is None:
        w = max(10, data.shape[1] * 0.35)
        h = max(6, data.shape[0] * 0.25)
        figsize = (min(w, 40), min(h, 60))

    g = sns.clustermap(
        data,
        cmap=cmap,
        z_score=z_score,
        figsize=figsize,
        xticklabels=True,
        yticklabels=True,
        linewidths=0.3,
        linecolor="white",
        dendrogram_ratio=(0.08, 0.12),
        cbar_kws={"label": "z-score" if z_score is not None else "Expression"},
    )
    g.fig.suptitle(title, y=1.02, fontsize=13, fontweight="bold")
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=7)

    if save_path:
        g.savefig(save_path, bbox_inches="tight")
        logger.info("Heatmap saved to %s", save_path)

    return g.fig


# ---------------------------------------------------------------------------
# Dot plot (expression level + detection rate)
# ---------------------------------------------------------------------------

def plot_dot_plot(
    mean_expr: pd.DataFrame,
    det_rate: pd.DataFrame,
    genes: list[str],
    cell_types: Optional[list[str]] = None,
    title: str = "Expression Dot Plot",
    figsize: Optional[tuple[float, float]] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Dot plot: dot size = detection rate, colour = mean expression.

    Parameters
    ----------
    mean_expr : pd.DataFrame
        Mean CPM matrix.
    det_rate : pd.DataFrame
        Detection-rate matrix.
    genes : list[str]
        Genes to display (y-axis).
    cell_types : list[str] | None
        Cell types to display (x-axis).  *None* → all.
    """
    if cell_types is None:
        cell_types = list(mean_expr.columns)

    genes_present = [g for g in genes if g in mean_expr.index]
    expr_sub = mean_expr.loc[genes_present, cell_types]
    det_sub = det_rate.loc[genes_present, cell_types]

    if figsize is None:
        figsize = (max(8, len(cell_types) * 0.4), max(4, len(genes_present) * 0.35))

    fig, ax = plt.subplots(figsize=figsize)

    # Build grid
    for yi, gene in enumerate(genes_present):
        for xi, ct in enumerate(cell_types):
            size = det_sub.loc[gene, ct] * 300  # scale
            colour_val = np.log2(1 + expr_sub.loc[gene, ct])
            ax.scatter(
                xi, yi,
                s=max(size, 1),
                c=colour_val,
                cmap="Reds",
                vmin=0,
                vmax=np.log2(1 + expr_sub.values.max()),
                edgecolors="grey",
                linewidths=0.3,
            )

    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(genes_present)))
    ax.set_yticklabels(genes_present, fontsize=7)
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Gene")

    # Add size legend
    for frac in [0.1, 0.25, 0.5, 1.0]:
        ax.scatter([], [], s=frac * 300, c="grey", alpha=0.5,
                   label=f"{frac:.0%} det.")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False,
              title="Detection Rate", fontsize=7, title_fontsize=8)

    # Colorbar
    sm = plt.cm.ScalarMappable(
        cmap="Reds",
        norm=plt.Normalize(0, np.log2(1 + expr_sub.values.max())),
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.15)
    cbar.set_label("log2(1 + CPM)", fontsize=8)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        logger.info("Dot plot saved to %s", save_path)

    return fig


# ---------------------------------------------------------------------------
# Single-gene expression profile across cell types
# ---------------------------------------------------------------------------

def plot_gene_profile(
    gene: str,
    mean_expr: pd.DataFrame,
    det_rate: pd.DataFrame,
    top_n: int = 20,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Bar chart of one gene's expression across its top-expressing cell types.

    Parameters
    ----------
    gene : str
        Gene symbol.
    mean_expr : pd.DataFrame
        Mean CPM matrix.
    det_rate : pd.DataFrame
        Detection-rate matrix.
    top_n : int
        Show only the top *N* cell types by expression.
    """
    if gene not in mean_expr.index:
        raise ValueError(f"Gene '{gene}' not found in expression matrix")

    expr = mean_expr.loc[gene].sort_values(ascending=False).head(top_n)
    det = det_rate.loc[gene].reindex(expr.index)

    fig, ax1 = plt.subplots(figsize=(10, 5))

    x = range(len(expr))
    bars = ax1.bar(x, expr.values, color="steelblue", alpha=0.8, edgecolor="white")
    ax1.set_xticks(x)
    ax1.set_xticklabels(expr.index, rotation=45, ha="right", fontsize=7)
    ax1.set_ylabel("Mean CPM", color="steelblue")
    ax1.set_title(title or f"{gene} — Expression Profile", fontweight="bold")

    # Overlay detection rate
    ax2 = ax1.twinx()
    ax2.plot(x, det.values, "ro-", markersize=4, linewidth=1, alpha=0.7)
    ax2.set_ylabel("Detection Rate", color="red")
    ax2.set_ylim(0, 1.05)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        logger.info("Gene profile saved to %s", save_path)

    return fig


# ---------------------------------------------------------------------------
# Specificity overview: Tau distribution
# ---------------------------------------------------------------------------

def plot_tau_distribution(
    tau: pd.Series,
    threshold: float = 0.85,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Histogram of Tau specificity index values.

    Parameters
    ----------
    tau : pd.Series
        Tau values per gene.
    threshold : float
        Vertical line at the specificity threshold.
    """
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(tau.dropna(), bins=100, color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(threshold, color="red", linestyle="--", linewidth=1.5,
               label=f"Threshold = {threshold}")
    n_pass = (tau >= threshold).sum()
    ax.set_xlabel("Tau Specificity Index")
    ax.set_ylabel("Number of Genes")
    ax.set_title(f"Tau Distribution — {n_pass:,} genes ≥ {threshold}", fontweight="bold")
    ax.legend()

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# Volcano-style: fold-change vs. detection rate
# ---------------------------------------------------------------------------

def plot_specificity_volcano(
    results: pd.DataFrame,
    cell_type: str,
    fc_threshold: float = 10.0,
    det_threshold: float = 0.25,
    label_top_n: int = 10,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Scatter plot of fold-change vs. on-target detection for one cell type.

    Parameters
    ----------
    results : pd.DataFrame
        Output of ``specificity.score_specificity``.
    cell_type : str
        Cell type to plot.
    label_top_n : int
        Number of top genes to label.
    """
    sub = results[results["cell_type"] == cell_type].copy()
    if sub.empty:
        logger.warning("No results for cell type '%s'", cell_type)
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, f"No results for {cell_type}", ha="center", va="center")
        return fig

    sub["log2_fc"] = np.log2(sub["fold_change"].clip(lower=0.01))

    fig, ax = plt.subplots(figsize=(8, 6))

    colours = sub["all_pass"].map({True: "red", False: "grey"})
    ax.scatter(sub["log2_fc"], sub["detection_rate"],
               c=colours, s=15, alpha=0.6, edgecolors="none")

    ax.axvline(np.log2(fc_threshold), color="blue", linestyle="--", alpha=0.5,
               label=f"FC = {fc_threshold}x")
    ax.axhline(det_threshold, color="green", linestyle="--", alpha=0.5,
               label=f"Det = {det_threshold}")

    # Label top hits
    top = sub[sub["all_pass"]].head(label_top_n)
    for _, row in top.iterrows():
        ax.annotate(
            row["gene"],
            (row["log2_fc"], row["detection_rate"]),
            fontsize=6, alpha=0.8,
            xytext=(5, 5), textcoords="offset points",
        )

    ax.set_xlabel("log₂(Fold-Change vs. Next Highest)")
    ax.set_ylabel("On-Target Detection Rate")
    ax.set_title(f"{cell_type} — Specificity Volcano", fontweight="bold")
    ax.legend(fontsize=8)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig
