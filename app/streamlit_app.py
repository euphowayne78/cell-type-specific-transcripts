"""Streamlit dashboard for querying cell-type specific transcripts."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.query import (
    get_available_cell_types,
    query_transcripts,
    compare_across_cell_types,
    generate_report,
)

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="Cell-Type Specific Transcripts",
    page_icon="🧬",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Data loading (cached)
# ---------------------------------------------------------------------------
DATA_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"


@st.cache_data
def load_data():
    """Load all required data files."""
    results_path = RESULTS_DIR / "specificity_scores.csv"
    mean_path = DATA_DIR / "preprocessed_mean_cpm.parquet"
    det_path = DATA_DIR / "preprocessed_detection_rate.parquet"

    missing = []
    if not results_path.exists():
        missing.append(str(results_path))
    if not mean_path.exists():
        missing.append(str(mean_path))
    if not det_path.exists():
        missing.append(str(det_path))

    if missing:
        return None, None, None, missing

    results = pd.read_csv(results_path)
    mean_df = pd.read_parquet(mean_path)
    det_df = pd.read_parquet(det_path)
    return results, mean_df, det_df, []


# ---------------------------------------------------------------------------
# Main app
# ---------------------------------------------------------------------------

def main():
    st.title("Cell-Type Specific Transcript Finder")
    st.markdown(
        "Identify transcripts uniquely expressed in specific human cell types — "
        "for use as ADAR sensor triggers in cell-type-restricted therapeutic mRNA."
    )

    # Load data
    results, mean_df, det_df, missing = load_data()

    if missing:
        st.error(
            "Required data files not found. Run the analysis notebooks first:\n\n"
            + "\n".join(f"- `{m}`" for m in missing)
        )
        st.info(
            "Run notebooks 01-03 in order:\n"
            "1. `01_data_acquisition.ipynb` — Download data from CELLxGENE Census\n"
            "2. `02_preprocessing.ipynb` — Filter and normalize\n"
            "3. `03_specificity_analysis.ipynb` — Compute specificity scores"
        )
        return

    # -----------------------------------------------------------------------
    # Sidebar: query controls
    # -----------------------------------------------------------------------
    st.sidebar.header("Query Parameters")

    available_cts = get_available_cell_types(results)
    selected_ct = st.sidebar.selectbox(
        "Select Cell Type",
        options=available_cts,
        index=0,
    )

    st.sidebar.markdown("---")
    st.sidebar.subheader("Threshold Controls")

    tau_min = st.sidebar.slider("Min Tau", 0.50, 1.00, 0.85, 0.01)
    fc_min = st.sidebar.slider("Min Fold-Change", 1.0, 100.0, 10.0, 1.0)
    det_min = st.sidebar.slider("Min Detection Rate", 0.0, 1.0, 0.25, 0.05)
    off_target_max = st.sidebar.slider("Max Off-Target Detection", 0.0, 0.50, 0.05, 0.01)
    expr_min = st.sidebar.slider("Min Mean CPM", 0.0, 50.0, 1.0, 0.5)
    top_n = st.sidebar.slider("Max Results", 10, 200, 50, 10)

    # -----------------------------------------------------------------------
    # Query
    # -----------------------------------------------------------------------
    hits = query_transcripts(
        results,
        cell_type=selected_ct,
        tau_min=tau_min,
        fc_min=fc_min,
        det_min=det_min,
        off_target_max=off_target_max,
        expr_min=expr_min,
        top_n=top_n,
    )

    # -----------------------------------------------------------------------
    # Results
    # -----------------------------------------------------------------------
    st.header(f"Results for: {selected_ct}")

    col1, col2, col3 = st.columns(3)
    col1.metric("Specific Transcripts", len(hits))
    col2.metric("Best Tau", f"{hits['tau'].max():.3f}" if not hits.empty else "N/A")
    col3.metric("Best Fold-Change", f"{hits['fold_change'].max():.1f}x" if not hits.empty else "N/A")

    if hits.empty:
        st.warning(
            f"No transcripts found for **{selected_ct}** with current thresholds. "
            "Try relaxing the threshold controls in the sidebar."
        )
        return

    # Results table
    st.subheader("Ranked Transcripts")
    display_cols = [
        "gene", "tau", "mean_cpm", "detection_rate",
        "fold_change", "max_off_target_det",
    ]
    if "annotation_tier" in hits.columns:
        display_cols.append("annotation_tier")

    st.dataframe(
        hits[display_cols].style.format({
            "tau": "{:.3f}",
            "mean_cpm": "{:.1f}",
            "detection_rate": "{:.1%}",
            "fold_change": "{:.1f}",
            "max_off_target_det": "{:.3f}",
        }),
        use_container_width=True,
        height=400,
    )

    # Download button
    csv = hits.to_csv(index=False)
    st.download_button(
        label="Download Results (CSV)",
        data=csv,
        file_name=f"specific_transcripts_{selected_ct.replace(' ', '_')}.csv",
        mime="text/csv",
    )

    # -----------------------------------------------------------------------
    # Visualizations
    # -----------------------------------------------------------------------
    st.markdown("---")
    st.subheader("Expression Heatmap")

    heatmap_genes = hits["gene"].head(min(30, len(hits))).tolist()
    heatmap_cts = list(mean_df.columns)

    # Build heatmap data
    heatmap_data = np.log2(1 + mean_df.loc[
        mean_df.index.isin(heatmap_genes), :
    ])

    if not heatmap_data.empty:
        # Z-score per gene
        z_data = heatmap_data.subtract(heatmap_data.mean(axis=1), axis=0)
        z_std = heatmap_data.std(axis=1)
        z_std[z_std == 0] = 1
        z_data = z_data.div(z_std, axis=0)

        fig = px.imshow(
            z_data.values,
            x=list(z_data.columns),
            y=list(z_data.index),
            color_continuous_scale="RdYlBu_r",
            labels=dict(x="Cell Type", y="Gene", color="Z-score"),
            aspect="auto",
        )
        fig.update_layout(
            height=max(300, len(heatmap_genes) * 22),
            xaxis_tickangle=45,
            font_size=9,
        )
        st.plotly_chart(fig, use_container_width=True)

    # -----------------------------------------------------------------------
    # Gene drill-down
    # -----------------------------------------------------------------------
    st.markdown("---")
    st.subheader("Gene Expression Profile")

    selected_gene = st.selectbox(
        "Select a gene to view its expression across all cell types:",
        options=hits["gene"].tolist(),
    )

    if selected_gene:
        comparison = compare_across_cell_types(results, selected_gene, mean_df, det_df)

        col_a, col_b = st.columns(2)

        with col_a:
            fig_bar = px.bar(
                comparison.head(25),
                x="cell_type",
                y="mean_cpm",
                title=f"{selected_gene} — Mean CPM (Top 25 Cell Types)",
                labels={"mean_cpm": "Mean CPM", "cell_type": "Cell Type"},
            )
            fig_bar.update_layout(xaxis_tickangle=45, height=400)
            st.plotly_chart(fig_bar, use_container_width=True)

        with col_b:
            fig_det = px.bar(
                comparison.head(25),
                x="cell_type",
                y="detection_rate",
                title=f"{selected_gene} — Detection Rate (Top 25)",
                labels={"detection_rate": "Detection Rate", "cell_type": "Cell Type"},
                color_discrete_sequence=["coral"],
            )
            fig_det.update_layout(xaxis_tickangle=45, height=400)
            st.plotly_chart(fig_det, use_container_width=True)

    # -----------------------------------------------------------------------
    # Text report
    # -----------------------------------------------------------------------
    st.markdown("---")
    with st.expander("View Text Report"):
        report = generate_report(results, selected_ct, mean_df, det_df, top_n=20)
        st.code(report, language=None)


if __name__ == "__main__":
    main()
