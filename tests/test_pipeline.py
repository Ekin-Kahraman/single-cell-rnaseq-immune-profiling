"""Tests for the single-cell RNA-seq pipeline."""

import scanpy as sc
import numpy as np
from pathlib import Path

# Import pipeline functions by path since filenames start with numbers
import importlib.util


def _import_script(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent.parent / "scripts" / f"{name}.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


qc_mod = _import_script("01_load_and_qc")
preprocess_mod = _import_script("02_preprocess")
cluster_mod = _import_script("04_cluster")
annotate_mod = _import_script("05_annotate_cell_types")


class TestQC:
    def test_qc_filters_cells(self, pbmc3k_raw, tmp_path):
        """QC should remove low-quality cells."""
        adata = pbmc3k_raw.copy()
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        n_before = adata.n_obs
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata = adata[adata.obs["n_genes_by_counts"] < 2500, :].copy()
        adata = adata[adata.obs["pct_counts_mt"] < 5, :].copy()

        assert adata.n_obs < n_before, "QC should remove some cells"
        assert adata.n_obs > 0, "QC should not remove all cells"

    def test_mito_annotation(self, pbmc3k_raw):
        """MT genes should be detected."""
        adata = pbmc3k_raw.copy()
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        assert adata.var["mt"].sum() > 0, "Should find mitochondrial genes"


class TestPreprocessing:
    def test_normalization(self, pbmc3k_raw):
        """Normalization should produce ~10k counts per cell."""
        adata = pbmc3k_raw.copy()
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)

        totals = np.array(adata.X.sum(axis=1)).flatten()
        np.testing.assert_allclose(totals, 1e4, rtol=1e-5)

    def test_hvg_selection(self, pbmc3k_raw):
        """HVG selection should identify a subset of genes."""
        adata = pbmc3k_raw.copy()
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)

        n_hvg = adata.var["highly_variable"].sum()
        assert 1000 <= n_hvg <= 2500, f"Expected ~2000 HVGs, got {n_hvg}"


class TestClustering:
    def test_leiden_produces_clusters(self, pbmc3k_raw):
        """Leiden should produce multiple clusters."""
        adata = pbmc3k_raw.copy()
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        adata = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=40)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
        sc.tl.leiden(adata, resolution=0.5)

        n_clusters = adata.obs["leiden"].nunique()
        assert n_clusters >= 3, f"Expected at least 3 clusters, got {n_clusters}"
        assert n_clusters <= 20, f"Too many clusters: {n_clusters}"


class TestAnnotation:
    def test_marker_dict_has_known_types(self):
        """PBMC marker dictionary should contain expected cell types."""
        expected = ["CD4+ T cells", "CD8+ T cells", "NK cells", "B cells", "CD14+ Monocytes"]
        for ct in expected:
            assert ct in annotate_mod.PBMC_MARKERS, f"Missing {ct}"

    def test_markers_are_valid_genes(self, pbmc3k_raw):
        """At least some marker genes should be present in the dataset."""
        all_markers = set()
        for genes in annotate_mod.PBMC_MARKERS.values():
            all_markers.update(genes)
        present = all_markers & set(pbmc3k_raw.var_names)
        assert len(present) >= 10, f"Only {len(present)} marker genes found in dataset"
