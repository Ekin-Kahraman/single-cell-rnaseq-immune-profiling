# Single-Cell RNA-seq Immune Cell Profiling

[![CI](https://github.com/Ekin-Kahraman/single-cell-rnaseq-immune-profiling/actions/workflows/ci.yml/badge.svg)](https://github.com/Ekin-Kahraman/single-cell-rnaseq-immune-profiling/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)

Which immune-cell populations can be identified from gene expression in a blood
sample? This reproducible [Scanpy](https://scanpy.readthedocs.io/) workflow takes
single-cell RNA data through quality control, clustering, marker-based cell
identification and closer analysis of T-cell subtypes.

The recorded run retains **2,604 cells**, assigns **five broad immune populations**,
and separates the T-cell cluster into candidate CD4/CD8 subtypes. It uses one
healthy donor and does not establish population-level or disease findings.

## Validation and reproducibility

- Automated tests run on Python 3.10, 3.11 and 3.12.
- A full pipeline check validates the cell counts, annotations, T-cell labels, result tables and figures.
- Each run records file checksums so outputs can be traced and checked. These software checks do not independently validate the biological labels.

<p align="center">
  <img src="docs/umap_3d_rotation.gif" alt="3D UMAP rotation showing PBMC immune cell clusters" width="600">
</p>

<details>
<summary>Analysis summary figure</summary>

![Publication Figure](docs/publication_figure.png)

**Panel A** - UMAP coloured by unsupervised Leiden clusters. **Panel B** - Same embedding coloured by assigned cell type. **Panel C** - Cell type proportions. **Panel D** - Z-scored expression of canonical marker genes per cell type (red = high, blue = low). **Panel E** - Summary statistics.
</details>

## Dataset

**10X Genomics PBMC 3k** - 2,700 peripheral blood mononuclear cells from a healthy donor, sequenced on the Chromium platform. This is the standard benchmark dataset used across [scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html), [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html), and other single-cell frameworks.

- **Direct download**: [filtered gene-barcode matrices](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) (5.9 MB)
- **Reference**: Zheng et al. (2017) [Massively parallel digital transcriptional profiling of single cells](https://doi.org/10.1038/ncomms14049). *Nature Communications* 8, 14049.

## Workflow

```
PBMC 3k (10X Genomics, 2,700 cells)
    │
    ▼
 01 QC ──────────── Scrublet doublet detection → filter: 200 < genes < 2500, mito < 5%
    │                (36 doublets detected, 34 removed)
    ▼
 02 Preprocess ──── Normalise (10k), log1p, 2000 HVGs, regress, scale
    │
    ▼
 03 Reduce ──────── PCA (40 PCs) → kNN graph → UMAP (random_state=42)
    │
    ▼
 04 Cluster ─────── Leiden at 5 resolutions → silhouette selection (≥5 clusters)
    │
    ▼
 05 Annotate ────── Wilcoxon DE → score against PBMC marker signatures → 5 cell types
    │
    ▼
 06 Trajectory ──── PAGA graph abstraction → diffusion pseudotime (rooted in CD14+ mono)
    │
    ▼
 07 Subcluster ──── T cell compartment → resolve CD4+ (47.5%) and CD8+ (52.5%)
    │
    ▼
 08 Figures ─────── Multi-panel publication figure (PNG 300 DPI + PDF vector)
    │
    ▼
 09 Manifest ─────── SHA-256 manifest for generated analysis artefacts
```

## Pipeline

| Step | Script | What it does |
|------|--------|--------------|
| 01 | `01_load_and_qc.py` | Download PBMC 3k, Scrublet doublet detection, QC metrics, filter low-quality cells and doublets |
| 02 | `02_preprocess.py` | Normalise to 10k counts/cell, log-transform, select 2,000 HVGs, regress out confounders, scale |
| 03 | `03_reduce_dimensions.py` | PCA (40 components), k-nearest neighbour graph, UMAP embedding |
| 04 | `04_cluster.py` | Leiden clustering at 5 resolutions (0.3–1.2), silhouette evaluation, select best with ≥5 cluster floor |
| 05 | `05_annotate_cell_types.py` | Wilcoxon rank-sum DE, score clusters against curated PBMC signatures, assign cell types |
| 06 | `06_trajectory.py` | PAGA partition-based graph abstraction, PAGA-initialised UMAP, diffusion pseudotime |
| 07 | `07_t_cell_subclustering.py` | Extract T cell compartment, subcluster, resolve CD4+/CD8+ via marker scoring |
| 08 | `08_publication_figures.py` | Multi-panel figure with UMAP, composition, marker heatmap, summary (PNG + PDF) |
| 09 | `09_output_manifest.py` | Generate checksums and file-size manifest for analysis artefacts |

## Project Structure

```
single-cell-rnaseq-immune-profiling/
├── scripts/
│   ├── 01_load_and_qc.py          QC + Scrublet doublet detection
│   ├── 02_preprocess.py            Normalise, HVG, regress, scale
│   ├── 03_reduce_dimensions.py     PCA, kNN graph, UMAP
│   ├── 04_cluster.py               Multi-resolution Leiden + silhouette
│   ├── 05_annotate_cell_types.py   Marker DE + automated annotation
│   ├── 06_trajectory.py            PAGA + diffusion pseudotime
│   ├── 07_t_cell_subclustering.py  CD4+/CD8+ resolution
│   ├── 08_publication_figures.py   Multi-panel figure (PNG + PDF)
│   ├── 09_output_manifest.py        Output checksums and file sizes
│   ├── validate_outputs.py          Full-run output validator
│   └── palette.py                  Shared Okabe-Ito colourblind palette
├── tests/
│   └── test_pipeline.py            7 tests (QC, normalisation, clustering, markers)
├── docs/
│   ├── umap_3d_rotation.gif        Animated 3D UMAP
│   └── publication_figure.png      Static multi-panel figure
├── run_pipeline.py                 CLI runner with --from resume flag
├── pyproject.toml                  Dependencies + metadata
├── requirements-lock.txt           Pinned versions for reproducibility
└── CITATION.cff
```

Each script reads the previous step's `.h5ad` output from `results/`.

## Results

### Cell Type Composition

| Cell Type | Cells | % | Key Markers |
|-----------|-------|---|-------------|
| T cells (initial coarse CD4+ label) | 1,192 | 45.8 | CD3D, IL7R |
| CD14+ Monocytes | 636 | 24.4 | CD14, LYZ |
| NK cells | 410 | 15.7 | NKG7, GNLY |
| B cells | 330 | 12.7 | MS4A1, CD79A |
| Dendritic cells | 36 | 1.4 | FCER1A, CST3 |

2,604 cells retained after QC and doublet removal (from 2,700 raw). Clustering selected resolution 0.5 (5 clusters, silhouette 0.204).

### T Cell Subclustering

The initial marker-based annotation labels the 1,192-cell T-cell cluster as
CD4+. A separate analysis of that cluster gives the following candidate labels:

| Subtype | Cells | % of T cells |
|---------|-------|---|
| CD8+ T | 626 | 52.5 |
| CD4+ T | 566 | 47.5 |

Labels use CD8A/CD8B/GZMK/GZMA scores compared with IL7R/CD4/TCF7/LEF1.
The split is a result of this marker-based procedure, not an independently
measured CD4/CD8 ratio.

### Trajectory Inference

PAGA summarises similarities between cell groups. Diffusion pseudotime orders
cells relative to a chosen CD14+ monocyte as an exploratory starting point.
In this single cross-sectional sample, neither result establishes lineage,
developmental direction, or a monocyte-to-dendritic-cell transition.

### Biological Interpretation

The broad labels recover familiar blood immune populations, with T cells forming
the largest cluster. A small candidate dendritic-cell group contains 36 cells.
Labels depend on marker choice and clustering resolution; they have not been
confirmed by an independent assay. FCGR3A+ monocytes are not separately resolved
at the selected resolution.

Silhouette scores in single-cell data are typically low due to continuous rather than discrete cell states; the metric is used here for relative comparison between resolutions, not as an absolute quality measure.

## Quick Start

```bash
git clone https://github.com/Ekin-Kahraman/single-cell-rnaseq-immune-profiling.git
cd single-cell-rnaseq-immune-profiling
pip install -e .                  # or: pip install -r requirements-lock.txt
python run_pipeline.py            # full pipeline (~38s)
python run_pipeline.py --from 6   # resume from trajectory step
python scripts/validate_outputs.py
```

Use `requirements-lock.txt` to pin dependency versions. Fixed seeds and recorded
outputs support reproducibility, but hardware and library differences can still
affect numerical results.

## Testing

```bash
pip install -e ".[dev]"
pytest -v
```

7 tests covering QC filtering, normalisation, HVG selection, clustering, and marker gene validation. CI runs tests on Python 3.10, 3.11, and 3.12, then runs and validates the full pipeline on Python 3.12.

## Design Decisions

- **Doublet detection** - Scrublet integrated before QC filtering. 36 doublets detected (1.3%), 34 removed after other QC filters. Recommended by [Luecken & Theis (2019)](https://doi.org/10.15252/msb.20188746).
- **Automated annotation** - Clusters scored against curated PBMC marker gene sets rather than manual inspection. The marker sets are themselves a subjective choice - but encoding them explicitly makes the annotation reproducible and auditable.
- **Multi-resolution clustering** - Leiden at 5 resolutions with silhouette evaluation. The ≥5 cluster floor reflects the known minimum of major PBMC lineages (T cells, B cells, monocytes, NK, DCs).
- **Exploratory cell connectivity** - PAGA summarises relationships between clusters. A CD14+ monocyte is selected as the pseudotime root for illustration; this does not make it a progenitor or establish a developmental pathway.
- **T cell subclustering** - Resolves CD4+/CD8+ populations that share CD3D/CD3E expression and cannot be separated at global clustering resolution.
- **Colourblind-friendly palette** - Okabe-Ito colours throughout.
- **Reproducible seeds** - `random_state=42` for UMAP, Leiden, Scrublet, and silhouette sampling.
- **Dual-format figures** - PNG (300 DPI) for web, PDF (vector) for publication submission.

## Limitations

- **Single-sample dataset.** Multi-sample work needs donor-aware analysis and assessment of technical differences; integration or batch correction should be justified for the study.
- **`regress_out` is debatable.** Used here following the original scanpy tutorial, but Luecken & Theis (2019) suggest regression may overcorrect for well-filtered cells.
- **No pathway enrichment.** Gene set enrichment (via decoupler or GSEApy) would connect cell types to functional programmes. Planned as a future addition.
- **FCGR3A+ monocytes not resolved.** Higher resolution or targeted analysis could be tested, but separation is not guaranteed.
- **Megakaryocytes not resolved.** The PBMC 3k dataset contains a small platelet/megakaryocyte population (PPBP+, PF4+) that merges with other clusters at this resolution. The canonical scanpy tutorial resolves 8 cell types from this dataset; our pipeline resolves 5 at the global level + 2 via T cell subclustering. The difference is resolution choice - we optimise for silhouette score rather than maximising cluster count.

## Licence

MIT
