"""Unit tests for specificity metrics."""

import numpy as np
import pandas as pd
import pytest

from src.specificity import (
    compute_tau,
    compute_fold_change,
    compute_detection_specificity,
    score_specificity,
    query_cell_type,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mock_expression():
    """Create a mock expression matrix with known specific genes."""
    np.random.seed(42)
    n_genes = 100
    n_types = 10
    cell_types = [f"cell_type_{i}" for i in range(n_types)]
    genes = [f"gene_{i}" for i in range(n_genes)]

    # Base: low random expression
    data = np.random.exponential(scale=0.5, size=(n_genes, n_types))

    # Plant specific genes:
    # gene_0: specific to cell_type_0 (high in ct0, low elsewhere)
    data[0, :] = 0.1
    data[0, 0] = 500.0

    # gene_1: specific to cell_type_1
    data[1, :] = 0.05
    data[1, 1] = 300.0

    # gene_2: ubiquitous (high in all)
    data[2, :] = 200.0

    # gene_3: not expressed
    data[3, :] = 0.0

    return pd.DataFrame(data, index=genes, columns=cell_types)


@pytest.fixture
def mock_detection():
    """Create a mock detection rate matrix matching mock_expression."""
    np.random.seed(42)
    n_genes = 100
    n_types = 10
    cell_types = [f"cell_type_{i}" for i in range(n_types)]
    genes = [f"gene_{i}" for i in range(n_genes)]

    data = np.random.uniform(0.01, 0.3, size=(n_genes, n_types))

    # gene_0: detected in ct0, not in others
    data[0, :] = 0.01
    data[0, 0] = 0.80

    # gene_1: detected in ct1, not in others
    data[1, :] = 0.02
    data[1, 1] = 0.60

    # gene_2: detected everywhere
    data[2, :] = 0.90

    # gene_3: not detected
    data[3, :] = 0.0

    return pd.DataFrame(data, index=genes, columns=cell_types)


# ---------------------------------------------------------------------------
# Tau tests
# ---------------------------------------------------------------------------

class TestComputeTau:
    def test_perfectly_specific_gene(self, mock_expression):
        tau = compute_tau(mock_expression)
        # gene_0 should have very high Tau (close to 1)
        assert tau.loc["gene_0"] > 0.95

    def test_ubiquitous_gene(self, mock_expression):
        tau = compute_tau(mock_expression)
        # gene_2 is expressed equally everywhere → Tau ≈ 0
        assert tau.loc["gene_2"] < 0.1

    def test_zero_expression_gene(self, mock_expression):
        tau = compute_tau(mock_expression)
        # gene_3 has zero expression → Tau = 0
        assert tau.loc["gene_3"] == 0.0

    def test_tau_range(self, mock_expression):
        tau = compute_tau(mock_expression)
        assert (tau >= 0).all()
        assert (tau <= 1).all()

    def test_output_is_series(self, mock_expression):
        tau = compute_tau(mock_expression)
        assert isinstance(tau, pd.Series)
        assert len(tau) == len(mock_expression)


# ---------------------------------------------------------------------------
# Fold-change tests
# ---------------------------------------------------------------------------

class TestFoldChange:
    def test_specific_gene_high_fc(self, mock_expression):
        fc_df, second_df = compute_fold_change(mock_expression)
        # gene_0 in cell_type_0 should have very high fold-change
        assert fc_df.loc["gene_0", "cell_type_0"] > 100

    def test_ubiquitous_gene_low_fc(self, mock_expression):
        fc_df, second_df = compute_fold_change(mock_expression)
        # gene_2 is ubiquitous → fold-change ≈ 1
        assert fc_df.loc["gene_2", "cell_type_0"] < 2

    def test_output_shape(self, mock_expression):
        fc_df, second_df = compute_fold_change(mock_expression)
        assert fc_df.shape == mock_expression.shape
        assert second_df.shape == mock_expression.shape


# ---------------------------------------------------------------------------
# Detection specificity tests
# ---------------------------------------------------------------------------

class TestDetectionSpecificity:
    def test_specific_gene_passes(self, mock_detection):
        det_pass = compute_detection_specificity(
            mock_detection, on_target_min=0.25, off_target_max=0.05,
        )
        # gene_0 in cell_type_0: 80% on-target, 1% off-target → True
        assert det_pass.loc["gene_0", "cell_type_0"] == True  # noqa: E712

    def test_ubiquitous_gene_fails(self, mock_detection):
        det_pass = compute_detection_specificity(
            mock_detection, on_target_min=0.25, off_target_max=0.05,
        )
        # gene_2 is detected everywhere → no cell type should pass
        assert not det_pass.loc["gene_2"].any()


# ---------------------------------------------------------------------------
# Full scoring pipeline
# ---------------------------------------------------------------------------

class TestScoreSpecificity:
    def test_recovers_planted_genes(self, mock_expression, mock_detection):
        results = score_specificity(
            mock_expression, mock_detection,
            tau_threshold=0.85,
            fc_threshold=10.0,
            on_target_det_min=0.25,
            off_target_det_max=0.05,
        )
        # Should find gene_0 and gene_1 as specific
        passing = results[results["all_pass"]]
        assert "gene_0" in passing["gene"].values
        assert "gene_1" in passing["gene"].values

    def test_ubiquitous_gene_excluded(self, mock_expression, mock_detection):
        results = score_specificity(
            mock_expression, mock_detection,
            tau_threshold=0.85,
        )
        # gene_2 (ubiquitous) should not appear (Tau too low)
        assert "gene_2" not in results["gene"].values

    def test_empty_with_extreme_thresholds(self, mock_expression, mock_detection):
        results = score_specificity(
            mock_expression, mock_detection,
            tau_threshold=0.9999,
            fc_threshold=1e6,
        )
        if not results.empty:
            assert results["all_pass"].sum() == 0

    def test_output_columns(self, mock_expression, mock_detection):
        results = score_specificity(mock_expression, mock_detection)
        if not results.empty:
            expected_cols = {
                "gene", "cell_type", "tau", "mean_cpm", "detection_rate",
                "fold_change", "all_pass",
            }
            assert expected_cols.issubset(set(results.columns))


# ---------------------------------------------------------------------------
# Query helper
# ---------------------------------------------------------------------------

class TestQueryCellType:
    def test_query_returns_correct_cell_type(self, mock_expression, mock_detection):
        results = score_specificity(mock_expression, mock_detection)
        if results.empty:
            pytest.skip("No results to query")
        hits = query_cell_type(results, "cell_type_0", all_pass_only=True)
        if not hits.empty:
            assert (hits["cell_type"] == "cell_type_0").all()

    def test_query_missing_cell_type(self, mock_expression, mock_detection):
        results = score_specificity(mock_expression, mock_detection)
        hits = query_cell_type(results, "nonexistent_cell_type")
        assert hits.empty
