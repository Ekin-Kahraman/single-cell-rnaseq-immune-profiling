"""Shared fixtures for pipeline tests."""

import pytest
import scanpy as sc


@pytest.fixture(scope="session")
def pbmc3k_raw():
    """Load the raw PBMC 3k dataset (cached across test session)."""
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    return adata
