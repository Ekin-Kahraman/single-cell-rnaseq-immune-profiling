# Reference run and verification

The public figures and small result tables were regenerated on **25 October 2026**
from the public 10x Genomics PBMC 3k dataset. The analysis used Python 3.12.13
on macOS arm64 with the exact packages in `requirements-lock.txt`.

## Reproduce

From the repository root, follow the README installation commands, then run:

```sh
python run_pipeline.py
python scripts/validate_outputs.py
pytest -v
ruff check run_pipeline.py scripts/ tests/
python scripts/make_3d_umap_gif.py
```

The full nine-stage analysis completed. Output validation, all **18 tests** and
lint passed. Tests cover both analysis primitives and the project's execution,
annotation and output-integrity contracts. The optional animation was generated
from the same annotated cells.

The subsequent display-label correction adds five regression tests (**23 total**).
The 3D view now joins the 638 CD4-like and 559 CD8-like labels by cell ID from
`07_t_cells.h5ad`; missing, duplicated or unrelated cells stop figure generation.
The static summary shows the six broad populations and labels the mixed
T-cell compartment simply as **T cells**. Neither display change alters the
underlying annotation scores, clusters or subtype counts.

The raw cached `.h5ad` input has SHA-256:

```text
f5ddd76da819fa241d7637754472cbf9067ebf18a27592383fec5df563d43a12
```

This hashes the Scanpy input object, not the original 10x archive.
Every completed run records its installed packages, platform, code revision,
working-tree state, executed stages and input checksum in
`results/run_metadata.json`. The output manifest hashes the generated analysis
files. The run record is finalised after that manifest and is not included in it.

## Reference results

| Measurement | Observed result |
| --- | ---: |
| Input cells | 2,700 |
| Retained cells | 2,604 |
| Highly variable genes | 2,000 |
| Selected Leiden resolution | 0.5 |
| Clusters / candidate broad populations | 6 |
| Silhouette score | 0.2188 |
| T-cell compartment | 1,197 |
| Candidate CD4-like / CD8-like cells | 638 / 559 |

[Machine-readable results](reference_results.json) and
[population counts](cell_type_composition.csv) accompany the updated
[summary figure](publication_figure.png).

## What changed from the previous presentation

The previous README reported five broad groups and 1,192 T cells, with a
566/626 CD4/CD8 split. A fresh run of the previous analysis code in the repaired,
pinned environment already produced the six-group result above. The workflow
fixes preserved those new reference results; they did not create a new biological
finding. The published figures and numbers now describe the environment that
can actually be installed and rerun.

The environment repair replaced an unavailable `seaborn==0.13.3` requirement
with 0.13.2 and resolved a complete Python 3.12 lockfile. The analysis now also
rejects missing annotation evidence, includes both broad T-cell labels in
subclustering, sizes PCA to the data, validates per-population counts and records
failed/resumed runs explicitly. Installing the repository still works after
`data/` and `results/` have been generated.

## Scope of verification

These results are marker-based candidates from one healthy donor. They are not
independently validated cell identities or disease findings. Fixed seeds and
pinned packages improve repeatability, but numerical libraries and platforms can
change clustering and labels; exact cell counts are a reference, not a cross-platform
test requirement.

[GitHub Actions](https://github.com/Ekin-Kahraman/single-cell-rnaseq-immune-profiling/actions/workflows/ci.yml)
tests Python 3.10–3.12 and runs the full dataset on Linux with the pinned Python
3.12 environment. Each full-run job uploads the run record, manifest, population
table and summary figure as an `analysis-evidence` artifact, keeping its own
observed results attached to the tested commit.
