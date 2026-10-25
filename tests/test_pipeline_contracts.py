"""Regression checks for project behaviour, separate from Scanpy's own tests."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_module(relative):
    spec = importlib.util.spec_from_file_location(Path(relative).stem, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = load_module("run_pipeline.py")
subtypes = load_module("scripts/07_t_cell_subclustering.py")
annotation = load_module("scripts/05_annotate_cell_types.py")
reduction = load_module("scripts/03_reduce_dimensions.py")
validator = load_module("scripts/validate_outputs.py")


@pytest.mark.parametrize("step", [0, -1, runner.N_STEPS + 1])
def test_invalid_resume_cannot_report_success(step, tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "PROJECT_DIR", tmp_path)
    with pytest.raises(ValueError, match="Start step"):
        runner.run_pipeline(step)
    assert not (tmp_path / "results").exists()


def test_failed_stage_stops_pipeline_and_records_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "PROJECT_DIR", tmp_path)
    monkeypatch.setattr(runner, "git_output", lambda *args: None)
    monkeypatch.setattr(runner.platform, "platform", lambda: "test platform")
    calls = []

    def execute(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(returncode=23)

    monkeypatch.setattr(runner.subprocess, "run", execute)
    with pytest.raises(SystemExit) as error:
        runner.run_pipeline(2)
    assert error.value.code == 23
    assert len(calls) == 1
    record = json.loads((tmp_path / "results/run_metadata.json").read_text())
    assert record["status"] == "failed"
    assert record["failed_step"] == 2
    assert record["completed_steps"] == []


def test_resume_record_identifies_the_steps_actually_executed(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "PROJECT_DIR", tmp_path)
    monkeypatch.setattr(runner, "git_output", lambda *args: None)
    monkeypatch.setattr(runner.platform, "platform", lambda: "test platform")
    monkeypatch.setattr(runner.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=0))
    runner.run_pipeline(runner.N_STEPS)
    record = json.loads((tmp_path / "results/run_metadata.json").read_text())
    assert record["status"] == "complete"
    assert record["start_step"] == runner.N_STEPS
    assert [step["step"] for step in record["completed_steps"]] == [runner.N_STEPS]


def t_cell_input():
    genes = ["IL7R", "CD4", "TCF7", "LEF1", "CD8A", "CD8B", "GZMK", "GZMA"]
    data = ad.AnnData(np.ones((10, len(genes))), var=pd.DataFrame(index=genes))
    data.obs["cell_type"] = ["CD4+ T cells"] * 4 + ["CD8+ T cells"] * 4 + ["B cells"] * 2
    data.raw = data.copy()
    data.obsm["X_pca"] = np.ones((10, 3))
    return data


def test_subclustering_keeps_cd8_as_well_as_cd4_parent_labels(monkeypatch):
    data = t_cell_input()
    monkeypatch.setattr(subtypes.sc.pp, "neighbors", lambda *args, **kwargs: None)
    monkeypatch.setattr(subtypes.sc.tl, "umap", lambda *args, **kwargs: None)
    monkeypatch.setattr(subtypes.sc.tl, "leiden", lambda *args, **kwargs: None)

    def score(adata, gene_list, score_name, **kwargs):
        adata.obs[score_name] = 1.0 if score_name == "cd4_score" else 0.0

    monkeypatch.setattr(subtypes.sc.tl, "score_genes", score)
    result = subtypes.subcluster_t_cells(data)
    assert result.obs_names.tolist() == data.obs_names[:8].tolist()
    assert set(result.obs["cell_type"]) == {"CD4+ T cells", "CD8+ T cells"}


def test_missing_t_cells_and_markers_fail_clearly():
    data = t_cell_input()
    data.obs["cell_type"] = "B cells"
    with pytest.raises(ValueError, match="four annotated T cells"):
        subtypes.subcluster_t_cells(data)
    data = t_cell_input()
    data.raw = data[:, ["CD4", "IL7R"]].copy()
    with pytest.raises(ValueError, match="Both CD4 and CD8"):
        subtypes.subcluster_t_cells(data)


def test_no_eligible_cell_type_signature_does_not_silently_assign_null_labels():
    data = ad.AnnData(np.ones((4, 3)), var=pd.DataFrame(index=["A", "B", "C"]))
    data.obs["leiden"] = pd.Categorical(["0", "0", "1", "1"])
    data.raw = data.copy()
    with pytest.raises(ValueError, match="No cell-type signature"):
        annotation.score_cell_types(data)


def test_small_feature_set_uses_available_pcs(tmp_path, monkeypatch):
    data = ad.AnnData(np.random.default_rng(42).normal(size=(20, 6)))
    monkeypatch.setattr(reduction, "FIG_DIR", tmp_path)
    # UMAP itself is exercised by the full real-data pipeline in CI.
    monkeypatch.setattr(reduction.sc.tl, "umap", lambda *args, **kwargs: None)
    result = reduction.reduce_dimensions(data)
    assert result.obsm["X_pca"].shape == (20, 5)
    assert result.uns["n_pcs_selected"] == 5
    assert result.uns["neighbors"]["params"]["n_pcs"] == 5


def test_composition_must_match_each_label_not_just_total(tmp_path, monkeypatch):
    data = ad.AnnData(np.ones((5, 2)))
    data.obs["cell_type"] = ["B"] * 3 + ["T"] * 2
    data.write_h5ad(tmp_path / "05_annotated.h5ad")
    pd.DataFrame({"count": [2, 3]}, index=["B", "T"]).to_csv(tmp_path / "cell_type_composition.csv")
    monkeypatch.setattr(validator, "RESULTS_DIR", tmp_path)
    with pytest.raises(SystemExit, match="each cell type"):
        validator.validate_tables()


def test_composition_matches_labels_regardless_of_category_order(tmp_path, monkeypatch):
    data = ad.AnnData(np.ones((5, 2)))
    data.obs["cell_type"] = pd.Categorical(["B"] * 3 + ["T"] * 2, categories=["T", "B"])
    data.write_h5ad(tmp_path / "05_annotated.h5ad")
    pd.DataFrame({"count": [3, 2]}, index=["B", "T"]).to_csv(tmp_path / "cell_type_composition.csv")
    pd.DataFrame({"group": ["B"] * 100, "names": ["CD79A"] * 100,
                  "scores": [1.0] * 100, "pvals_adj": [0.01] * 100,
                  "logfoldchanges": [2.0] * 100}).to_csv(tmp_path / "marker_genes.csv", index=False)
    monkeypatch.setattr(validator, "RESULTS_DIR", tmp_path)
    validator.validate_tables()
