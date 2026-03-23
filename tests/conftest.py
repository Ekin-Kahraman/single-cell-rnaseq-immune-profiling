"""Shared fixtures for pipeline tests."""

import pytest
import scanpy as sc
import numpy as np


@pytest.fixture(scope="session")
def pbmc3k_raw():
    """Load the raw PBMC 3k dataset (cached across test session)."""
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    return adata


@pytest.fixture
def small_adata():
    """Create a small synthetic AnnData for fast unit tests."""
    rng = np.random.default_rng(42)
    n_cells, n_genes = 200, 500
    X = rng.poisson(1, size=(n_cells, n_genes)).astype(np.float32)

    adata = sc.AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    # Add some "mitochondrial" genes
    mt_genes = [f"MT-{c}" for c in ["ND1", "ND2", "CO1", "CO2", "ATP6"]]
    for i, name in enumerate(mt_genes):
        if i < n_genes:
            adata.var_names = adata.var_names.tolist()
            adata.var_names = [name if j == i else adata.var_names[j] for j in range(n_genes)]

    adata.var_names_make_unique()
    return adata
