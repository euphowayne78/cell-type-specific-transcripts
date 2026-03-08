"""Unit tests for filtering modules."""

import numpy as np
import pandas as pd
import pytest

from src.preprocessing import (
    filter_low_expression,
    remove_housekeeping,
    flag_housekeeping_and_ribosomal,
    log1p_transform,
    run_preprocessing,
)
from src.filters import (
    annotate_quality_tiers,
    filter_cytoplasmic_lncrnas,
    filter_by_length,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def small_expression():
    """Small expression matrix for testing."""
    genes = ["ALB", "GAPDH", "RPL13A", "MT-CO1", "NOVEL1", "LOWEXPR"]
    cts = ["hepatocyte", "T cell", "neuron"]
    data = np.array([
        [500.0, 0.1, 0.05],   # ALB: hepatocyte-specific
        [200.0, 180.0, 190.0], # GAPDH: housekeeping
        [100.0, 90.0, 110.0],  # RPL13A: ribosomal
        [50.0, 40.0, 60.0],    # MT-CO1: mitochondrial
        [0.0, 300.0, 0.0],     # NOVEL1: T-cell specific
        [0.01, 0.005, 0.002],  # LOWEXPR: low expression
    ])
    return pd.DataFrame(data, index=genes, columns=cts)


@pytest.fixture
def small_detection():
    """Detection rate matrix matching small_expression."""
    genes = ["ALB", "GAPDH", "RPL13A", "MT-CO1", "NOVEL1", "LOWEXPR"]
    cts = ["hepatocyte", "T cell", "neuron"]
    data = np.array([
        [0.80, 0.01, 0.005],
        [0.95, 0.93, 0.94],
        [0.90, 0.88, 0.91],
        [0.70, 0.65, 0.72],
        [0.00, 0.50, 0.00],
        [0.005, 0.002, 0.001],
    ])
    return pd.DataFrame(data, index=genes, columns=cts)


@pytest.fixture
def mock_results():
    """Mock specificity results for filter tests."""
    return pd.DataFrame({
        "gene": ["ALB", "CD3D", "NOVEL1", "XIST", "SHORT1"],
        "cell_type": ["hepatocyte", "T cell", "T cell", "T cell", "neuron"],
        "tau": [0.95, 0.90, 0.88, 0.87, 0.86],
        "mean_cpm": [500, 200, 100, 50, 10],
        "detection_rate": [0.8, 0.6, 0.5, 0.3, 0.3],
        "fold_change": [100, 50, 30, 20, 15],
        "max_off_target_det": [0.01, 0.02, 0.03, 0.04, 0.04],
        "all_pass": [True, True, True, True, True],
    })


# ---------------------------------------------------------------------------
# Preprocessing tests
# ---------------------------------------------------------------------------

class TestFilterLowExpression:
    def test_removes_low_genes(self, small_expression, small_detection):
        mean_out, det_out = filter_low_expression(
            small_expression, small_detection, min_cpm=1.0, min_det_rate=0.01,
        )
        # LOWEXPR should be removed (max CPM = 0.01, below 1.0)
        assert "LOWEXPR" not in mean_out.index

    def test_keeps_expressed_genes(self, small_expression, small_detection):
        mean_out, det_out = filter_low_expression(
            small_expression, small_detection, min_cpm=1.0, min_det_rate=0.01,
        )
        assert "ALB" in mean_out.index
        assert "GAPDH" in mean_out.index


class TestRemoveHousekeeping:
    def test_removes_hk_genes(self, small_expression, small_detection):
        mean_out, det_out = remove_housekeeping(small_expression, small_detection)
        assert "GAPDH" not in mean_out.index

    def test_removes_ribosomal(self, small_expression, small_detection):
        mean_out, det_out = remove_housekeeping(small_expression, small_detection)
        assert "RPL13A" not in mean_out.index

    def test_removes_mitochondrial(self, small_expression, small_detection):
        mean_out, det_out = remove_housekeeping(small_expression, small_detection)
        assert "MT-CO1" not in mean_out.index

    def test_keeps_non_hk(self, small_expression, small_detection):
        mean_out, det_out = remove_housekeeping(small_expression, small_detection)
        assert "ALB" in mean_out.index
        assert "NOVEL1" in mean_out.index


class TestFlagHousekeeping:
    def test_flags_known_hk(self):
        genes = pd.Index(["GAPDH", "ACTB", "ALB", "CD3D"])
        flags = flag_housekeeping_and_ribosomal(genes)
        assert flags[0] == True   # GAPDH
        assert flags[1] == True   # ACTB
        assert flags[2] == False  # ALB
        assert flags[3] == False  # CD3D

    def test_flags_ribosomal(self):
        genes = pd.Index(["RPL13A", "RPS18", "MRPL1", "MRPS2", "ALB"])
        flags = flag_housekeeping_and_ribosomal(genes)
        assert flags[:4].all()
        assert flags[4] == False


class TestLog1pTransform:
    def test_transform_values(self, small_expression):
        result = log1p_transform(small_expression)
        # log2(1 + 500) ≈ 8.97
        assert abs(result.loc["ALB", "hepatocyte"] - np.log2(501)) < 0.01

    def test_zero_stays_zero(self):
        df = pd.DataFrame({"a": [0.0]}, index=["g1"])
        result = log1p_transform(df)
        assert result.iloc[0, 0] == 0.0


class TestRunPreprocessing:
    def test_pipeline(self, small_expression, small_detection):
        mean_out, det_out = run_preprocessing(
            small_expression, small_detection,
            min_cpm=1.0, min_det_rate=0.01,
        )
        # Should remove LOWEXPR, GAPDH, RPL13A, MT-CO1
        assert "LOWEXPR" not in mean_out.index
        assert "GAPDH" not in mean_out.index
        assert "ALB" in mean_out.index


# ---------------------------------------------------------------------------
# Post-filter tests
# ---------------------------------------------------------------------------

class TestAnnotateQualityTiers:
    def test_with_uniprot(self, mock_results):
        reviewed = {"ALB", "CD3D"}
        out = annotate_quality_tiers(mock_results, uniprot_reviewed=reviewed)
        assert out.loc[out["gene"] == "ALB", "annotation_tier"].iloc[0] == 1
        assert out.loc[out["gene"] == "NOVEL1", "annotation_tier"].iloc[0] == 2

    def test_without_uniprot(self, mock_results):
        out = annotate_quality_tiers(mock_results, uniprot_reviewed=None)
        assert (out["annotation_tier"] == 2).all()


class TestFilterCytoplasmicLncrnas:
    def test_removes_nuclear(self, mock_results):
        lncatlas = pd.DataFrame({
            "gene_name": ["XIST", "NOVEL1"],
            "RCI": [-2.0, 0.5],  # XIST nuclear, NOVEL1 cytoplasmic
        })
        out = filter_cytoplasmic_lncrnas(mock_results, lncatlas_df=lncatlas)
        assert "XIST" not in out["gene"].values
        assert "NOVEL1" in out["gene"].values

    def test_no_filter_without_data(self, mock_results):
        out = filter_cytoplasmic_lncrnas(mock_results, lncatlas_df=None)
        assert len(out) == len(mock_results)


class TestFilterByLength:
    def test_removes_short(self, mock_results):
        length_map = pd.Series(
            [2000, 1500, 800, 100, 50],
            index=["ALB", "CD3D", "NOVEL1", "XIST", "SHORT1"],
        )
        out = filter_by_length(mock_results, gene_length_map=length_map, min_length=200)
        assert "SHORT1" not in out["gene"].values
        assert "XIST" not in out["gene"].values
        assert "ALB" in out["gene"].values
