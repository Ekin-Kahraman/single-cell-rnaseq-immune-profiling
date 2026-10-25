"""Ensure prominent figures use the correct annotation level and cell identities."""

import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    "visual_labels", Path(__file__).resolve().parents[1] / "scripts/visual_labels.py"
)
labels = importlib.util.module_from_spec(spec)
spec.loader.exec_module(labels)


def example():
    parent = ad.AnnData(np.ones((4, 2)))
    parent.obs_names = ["a", "b", "c", "d"]
    parent.obs["cell_type"] = ["CD4+ T cells", "B cells", "CD4+ T cells", "CD8+ T cells"]
    subset = parent[["d", "a", "c"]].copy()
    subset.obs["t_subtype"] = ["CD8+ T", "CD4+ T", "CD8+ T"]
    return parent, subset


def test_refined_labels_follow_cell_ids_not_row_order():
    parent, subset = example()
    assert labels.refined_cell_types(parent, subset).tolist() == [
        "CD4-like T cells", "B cells", "CD8-like T cells", "CD8-like T cells",
    ]
    assert labels.broad_cell_types(parent).tolist() == ["T cells", "B cells", "T cells", "T cells"]
    assert parent.obs.loc["a", "cell_type"] == "CD4+ T cells"


@pytest.mark.parametrize("case", ["missing", "unrelated", "duplicate", "unknown_subtype"])
def test_incomplete_or_misaligned_refinement_cannot_generate_a_figure(case):
    parent, subset = example()
    if case == "missing":
        subset = subset[:2].copy()
    elif case == "unrelated":
        subset.obs_names = ["b", "a", "c"]
    elif case == "duplicate":
        subset.obs_names = ["a", "a", "c"]
    else:
        subset.obs.loc["a", "t_subtype"] = "unknown"
    with pytest.raises(ValueError):
        labels.refined_cell_types(parent, subset)
