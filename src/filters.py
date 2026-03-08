"""Post-specificity filtering: annotation quality, localization, transcript length."""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Annotation quality tiers
# ---------------------------------------------------------------------------

def annotate_quality_tiers(
    results: pd.DataFrame,
    uniprot_reviewed: Optional[set[str]] = None,
) -> pd.DataFrame:
    """Add an ``annotation_tier`` column to specificity results.

    Tiers:
    - **Tier 1**: Gene in UniProt reviewed set (Swiss-Prot) — well-characterised.
    - **Tier 2**: Gene not in reviewed set but has an ENSEMBL protein-coding
      annotation — predicted / unreviewed.
    - **Tier 3**: Everything else (novel / uncharacterised).

    Parameters
    ----------
    results : pd.DataFrame
        Output of ``specificity.score_specificity``.
    uniprot_reviewed : set[str] | None
        Set of gene symbols with UniProt reviewed (Swiss-Prot) entries.
        If *None*, all genes are assigned Tier 2.

    Returns
    -------
    pd.DataFrame
        Copy with ``annotation_tier`` column added.
    """
    out = results.copy()
    if uniprot_reviewed is None:
        out["annotation_tier"] = 2
        logger.info("No UniProt reviewed set provided — all genes assigned Tier 2")
        return out

    out["annotation_tier"] = out["gene"].apply(
        lambda g: 1 if g in uniprot_reviewed else 2
    )
    n_t1 = (out["annotation_tier"] == 1).sum()
    logger.info("Annotation tiers: Tier 1 = %d, Tier 2 = %d", n_t1, len(out) - n_t1)
    return out


def load_uniprot_reviewed_genes(filepath: Path | str) -> set[str]:
    """Load UniProt reviewed gene symbols from a TSV download.

    Expected file: a TSV with at least a ``Gene Names (primary)`` column,
    downloaded from https://www.uniprot.org/uniprot/?query=reviewed:true+AND+organism_id:9606

    Parameters
    ----------
    filepath : Path | str
        Path to the downloaded UniProt TSV.

    Returns
    -------
    set[str]
        Gene symbols.
    """
    df = pd.read_csv(filepath, sep="\t")
    col = None
    for candidate in ["Gene Names (primary)", "Gene names  (primary )", "Gene_Name"]:
        if candidate in df.columns:
            col = candidate
            break
    if col is None:
        raise ValueError(
            f"Could not find gene name column in {filepath}. "
            f"Available columns: {list(df.columns)}"
        )
    genes = set(df[col].dropna().str.strip())
    logger.info("Loaded %d UniProt reviewed gene symbols from %s", len(genes), filepath)
    return genes


# ---------------------------------------------------------------------------
# Subcellular localization filter (for lncRNAs)
# ---------------------------------------------------------------------------

def filter_cytoplasmic_lncrnas(
    results: pd.DataFrame,
    lncatlas_df: Optional[pd.DataFrame] = None,
    rci_threshold: float = 0.0,
) -> pd.DataFrame:
    """Remove lncRNAs that are not cytoplasmic.

    Uses the Relative Concentration Index (RCI) from lncATLAS:
    - RCI > 0 → cytoplasmic enrichment
    - RCI < 0 → nuclear enrichment

    Protein-coding genes are always kept (assumed cytosolic).

    Parameters
    ----------
    results : pd.DataFrame
        Specificity results.
    lncatlas_df : pd.DataFrame | None
        lncATLAS data with columns ``gene_name`` and ``RCI``.
        If *None*, no lncRNA filtering is performed.
    rci_threshold : float
        Minimum RCI to consider a lncRNA cytoplasmic.

    Returns
    -------
    pd.DataFrame
        Filtered results.
    """
    if lncatlas_df is None:
        logger.info("No lncATLAS data provided — skipping cytoplasmic filter")
        return results

    # Build RCI lookup
    rci_map = lncatlas_df.set_index("gene_name")["RCI"].to_dict()

    def is_cytoplasmic(gene: str) -> bool:
        if gene in rci_map:
            return rci_map[gene] >= rci_threshold
        # If gene not in lncATLAS, assume protein-coding → keep
        return True

    keep_mask = results["gene"].apply(is_cytoplasmic)
    n_removed = (~keep_mask).sum()
    logger.info("Cytoplasmic filter: removed %d nuclear lncRNAs (RCI < %.1f)", n_removed, rci_threshold)
    return results.loc[keep_mask].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Transcript length filter
# ---------------------------------------------------------------------------

def filter_by_length(
    results: pd.DataFrame,
    gene_length_map: Optional[pd.Series] = None,
    min_length: int = 200,
) -> pd.DataFrame:
    """Remove genes with transcripts shorter than *min_length* nt.

    Short transcripts are impractical for ADAR sensor design (need ~20-40 nt
    complementary region plus flanking context).

    Parameters
    ----------
    gene_length_map : pd.Series | None
        Series mapping gene name → transcript length.  If *None*, skips.
    min_length : int
        Minimum transcript length in nucleotides.
    """
    if gene_length_map is None:
        logger.info("No gene length data — skipping length filter")
        return results

    lengths = results["gene"].map(gene_length_map)
    keep = lengths.fillna(min_length) >= min_length  # keep unknowns
    n_removed = (~keep).sum()
    logger.info("Length filter (>= %d nt): removed %d genes", min_length, n_removed)
    return results.loc[keep].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Combined post-filtering pipeline
# ---------------------------------------------------------------------------

def apply_post_filters(
    results: pd.DataFrame,
    uniprot_reviewed: Optional[set[str]] = None,
    lncatlas_df: Optional[pd.DataFrame] = None,
    gene_length_map: Optional[pd.Series] = None,
    min_length: int = 200,
    rci_threshold: float = 0.0,
) -> pd.DataFrame:
    """Apply all post-specificity filters and annotation.

    1. Annotate quality tiers.
    2. Remove nuclear lncRNAs.
    3. Remove short transcripts.

    Returns
    -------
    pd.DataFrame
        Filtered and annotated results.
    """
    out = annotate_quality_tiers(results, uniprot_reviewed=uniprot_reviewed)
    out = filter_cytoplasmic_lncrnas(out, lncatlas_df=lncatlas_df, rci_threshold=rci_threshold)
    out = filter_by_length(out, gene_length_map=gene_length_map, min_length=min_length)
    logger.info("Post-filtering complete: %d genes remain", len(out))
    return out
