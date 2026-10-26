"""The checksum manifest must cover all required analysis artefacts once."""

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def manifest_outputs(tmp_path, monkeypatch):
    """Use small files to exercise integrity checks, not scientific analysis."""
    monkeypatch.chdir(tmp_path)
    validator = load_script("validate_outputs")
    writer = load_script("09_output_manifest")
    for name in validator.REQUIRED_FILES:
        if name in {"output_manifest.csv", "run_metadata.json"}:
            continue
        path = Path("results") / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"fixture for {name}\n".encode())
    writer.main()
    return validator, Path("results/output_manifest.csv")


def test_generated_manifest_covers_required_outputs(manifest_outputs):
    validator, _ = manifest_outputs
    validator.validate_manifest()


def test_manifest_cannot_omit_a_changed_required_output(manifest_outputs):
    validator, manifest_path = manifest_outputs
    omitted = "results/06_trajectory.h5ad"
    manifest = pd.read_csv(manifest_path)
    manifest.loc[manifest["path"] != omitted].to_csv(manifest_path, index=False)
    Path(omitted).write_bytes(b"changed after the manifest was created\n")
    with pytest.raises(SystemExit, match=r"omits required outputs.*06_trajectory.h5ad"):
        validator.validate_manifest()


@pytest.mark.parametrize("alias", [False, True], ids=["same-path", "equivalent-path"])
def test_manifest_rejects_duplicate_artefact_records(manifest_outputs, alias):
    validator, manifest_path = manifest_outputs
    manifest = pd.read_csv(manifest_path)
    duplicate = manifest.iloc[[0]].copy()
    if alias:
        duplicate["path"] = duplicate["path"].str.replace("results/", "results/./", regex=False)
    pd.concat([manifest, duplicate], ignore_index=True).to_csv(manifest_path, index=False)
    with pytest.raises(SystemExit, match="duplicate artefact paths"):
        validator.validate_manifest()
