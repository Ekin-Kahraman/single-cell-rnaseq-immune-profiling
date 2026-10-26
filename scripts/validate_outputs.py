#!/usr/bin/env python3
"""Validate outputs from a completed single-cell pipeline run."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import scanpy as sc

RESULTS_DIR = Path("results")

REQUIRED_FILES = [
    "01_filtered.h5ad",
    "02_preprocessed.h5ad",
    "03_reduced.h5ad",
    "04_clustered.h5ad",
    "05_annotated.h5ad",
    "06_trajectory.h5ad",
    "07_t_cells.h5ad",
    "marker_genes.csv",
    "cell_type_composition.csv",
    "figures/01_qc_metrics.png",
    "figures/02_hvg_selection.png",
    "figures/03_pca_variance.png",
    "figures/04_clustering.png",
    "figures/05_marker_dotplot.png",
    "figures/05_cell_types_umap.png",
    "figures/06_paga_graph.png",
    "figures/06_paga_graph.pdf",
    "figures/06_trajectory.png",
    "figures/06_trajectory.pdf",
    "figures/07_t_cell_subclustering.png",
    "figures/07_t_cell_subclustering.pdf",
    "figures/07_t_cell_markers.png",
    "figures/08_publication_figure.png",
    "figures/08_publication_figure.pdf",
    "output_manifest.csv",
    "run_metadata.json",
]


def require_file(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"Missing expected output: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"Output is empty: {path}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_manifest() -> None:
    manifest_path = RESULTS_DIR / "output_manifest.csv"
    manifest = pd.read_csv(manifest_path)
    required_cols = {"path", "bytes", "sha256"}
    if not required_cols.issubset(manifest.columns):
        raise SystemExit(f"Manifest missing columns: {required_cols - set(manifest.columns)}")

    paths = [Path(value).resolve() for value in manifest["path"]]
    if len(set(paths)) != len(paths):
        raise SystemExit("Manifest contains duplicate artefact paths")

    # The manifest cannot hash itself, and the runner finalises its metadata
    # after stage 09. Every other required analysis output must be covered.
    required_paths = {
        (RESULTS_DIR / name).resolve()
        for name in REQUIRED_FILES
        if name not in {"output_manifest.csv", "run_metadata.json"}
    }
    missing = sorted(required_paths - set(paths))
    if missing:
        raise SystemExit(f"Manifest omits required outputs: {', '.join(map(str, missing))}")

    for row in manifest.itertuples(index=False):
        path = Path(row.path)
        require_file(path)
        if int(row.bytes) != path.stat().st_size:
            raise SystemExit(f"Manifest size mismatch for {path}")
        if row.sha256 != sha256_file(path):
            raise SystemExit(f"Manifest checksum mismatch for {path}")


def validate_h5ad_outputs() -> None:
    adata = sc.read_h5ad(RESULTS_DIR / "06_trajectory.h5ad")
    if not 2500 <= adata.n_obs <= 2700:
        raise SystemExit(f"Unexpected retained cell count: {adata.n_obs}")
    for col in ["leiden", "cell_type", "dpt_pseudotime"]:
        if col not in adata.obs:
            raise SystemExit(f"Missing obs column in trajectory object: {col}")
    if adata.obs["cell_type"].nunique() < 5:
        raise SystemExit("Expected at least five annotated cell types")
    if adata.obs["cell_type"].isna().any():
        raise SystemExit("Some retained cells have no cell-type annotation")

    t_cells = sc.read_h5ad(RESULTS_DIR / "07_t_cells.h5ad")
    if t_cells.n_obs < 500:
        raise SystemExit(f"Unexpected T cell subset size: {t_cells.n_obs}")
    if "t_subtype" not in t_cells.obs:
        raise SystemExit("Missing T cell subtype annotations")


def validate_tables() -> None:
    composition = pd.read_csv(RESULTS_DIR / "cell_type_composition.csv", index_col=0)
    if composition.empty:
        raise SystemExit("Cell type composition table is empty")
    adata = sc.read_h5ad(RESULTS_DIR / "05_annotated.h5ad")
    expected = adata.obs["cell_type"].value_counts()
    expected.index = expected.index.astype(str)
    expected = expected.sort_index()
    reported = composition.iloc[:, 0].sort_index()
    if not reported.index.is_unique or not reported.index.equals(expected.index):
        raise SystemExit("Composition cell-type labels do not match annotated cells")
    if not (reported.to_numpy() == expected.to_numpy()).all():
        raise SystemExit("Composition counts do not match the annotations for each cell type")

    markers = pd.read_csv(RESULTS_DIR / "marker_genes.csv")
    required_marker_cols = {"group", "names", "scores", "pvals_adj", "logfoldchanges"}
    if not required_marker_cols.issubset(markers.columns):
        raise SystemExit("Marker table is missing expected rank_genes_groups columns")
    if len(markers) < 100:
        raise SystemExit("Marker table contains too few marker rows")


def main() -> None:
    for relpath in REQUIRED_FILES:
        require_file(RESULTS_DIR / relpath)
    validate_manifest()
    validate_h5ad_outputs()
    validate_tables()
    record = json.loads((RESULTS_DIR / "run_metadata.json").read_text())
    if record.get("status") != "complete":
        raise SystemExit("The most recent pipeline run did not complete successfully")
    print("Validated single-cell pipeline outputs")


if __name__ == "__main__":
    main()
