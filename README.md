# Single-Cell RNA-seq Immune Cell Profiling

[![CI](https://github.com/Ekin-Kahraman/single-cell-rnaseq-immune-profiling/actions/workflows/ci.yml/badge.svg)](https://github.com/Ekin-Kahraman/single-cell-rnaseq-immune-profiling/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)

End-to-end single-cell RNA-seq analysis pipeline in Python using [scanpy](https://scanpy.readthedocs.io/). Quality control, normalisation, dimensionality reduction, unsupervised clustering with automated resolution selection, and marker-based cell type annotation on human peripheral blood mononuclear cells.

<p align="center">
  <img src="docs/umap_3d_rotation.gif" alt="3D UMAP rotation showing PBMC immune cell clusters" width="600">
</p>

<details>
<summary>Static publication figure</summary>

![Publication Figure](docs/publication_figure.png)

**Panel A** — UMAP coloured by unsupervised Leiden clusters. **Panel B** — Same embedding coloured by assigned cell type. **Panel C** — Cell type proportions. **Panel D** — Z-scored expression of canonical marker genes per cell type (red = high, blue = low). **Panel E** — Summary statistics.
</details>

## Dataset

**10X Genomics PBMC 3k** — 2,700 peripheral blood mononuclear cells from a healthy donor, sequenced on the Chromium platform. This is the standard benchmark dataset used across [scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html), [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html), and other single-cell frameworks.

- **Direct download**: [filtered gene-barcode matrices](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) (5.9 MB)
- **Reference**: Zheng et al. (2017) [Massively parallel digital transcriptional profiling of single cells](https://doi.org/10.1038/ncomms14049). *Nature Communications* 8, 14049.

## Workflow

```
PBMC 3k (10X Genomics)
    │
    ▼
 01 QC ──────────── Filter: 200 < genes < 2500, mito < 5%
    │
    ▼
 02 Preprocess ──── Normalise (10k), log1p, 2000 HVGs, regress, scale
    │
    ▼
 03 Reduce ──────── PCA (40 PCs) → kNN graph → UMAP
    │
    ▼
 04 Cluster ─────── Leiden at 5 resolutions → silhouette selection (≥5 clusters)
    │
    ▼
 05 Annotate ────── Wilcoxon DE → score against PBMC marker signatures
    │
    ▼
 06 Figures ─────── Multi-panel publication figure + 3D UMAP
```

## Pipeline

| Step | Script | What it does |
|------|--------|--------------|
| 01 | `01_load_and_qc.py` | Download PBMC 3k, calculate QC metrics (genes/cell, UMI counts, mitochondrial %), filter low-quality cells |
| 02 | `02_preprocess.py` | Normalise to 10k counts/cell, log-transform, select 2,000 highly variable genes, regress out confounders, scale |
| 03 | `03_reduce_dimensions.py` | PCA (40 components), build k-nearest neighbour graph, compute UMAP embedding |
| 04 | `04_cluster.py` | Leiden clustering at 5 resolutions (0.3–1.2), evaluate with silhouette score, select best with a floor of 5 clusters |
| 05 | `05_annotate_cell_types.py` | Wilcoxon rank-sum test for marker genes, score clusters against known PBMC signatures, assign cell types |
| 06 | `06_publication_figures.py` | Multi-panel figure: UMAP, composition bar chart, marker heatmap, summary statistics |

All scripts are in `scripts/`. Each reads the previous step's `.h5ad` output from `results/`.

## Results

| Cell Type | Cells | % | Key Markers |
|-----------|-------|---|-------------|
| CD4+ T cells | 1,195 | 45.3 | CD3D, IL7R |
| CD14+ Monocytes | 464 | 17.6 | CD14, LYZ |
| NK cells | 419 | 15.9 | NKG7, GNLY |
| B cells | 342 | 13.0 | MS4A1, CD79A |
| FCGR3A+ Monocytes | 180 | 6.8 | FCGR3A, MS4A7 |
| Dendritic cells | 38 | 1.4 | FCER1A, CST3 |

The dominance of CD4+ T cells (45%) is expected in healthy donor PBMCs. The ratio of classical (CD14+) to nonclassical (FCGR3A+) monocytes is approximately 2.6:1, consistent with published literature. Dendritic cells are a rare population (1.4%), correctly resolved as a distinct cluster. CD8+ T cells and megakaryocytes are present in the dataset but were not resolved as separate clusters at resolution 0.5 — they likely merge with the CD4+ T cell and monocyte clusters respectively due to shared marker expression (CD3D/CD3E for T cell subtypes).

Clustering selected resolution 0.5 (6 clusters, silhouette 0.196). Silhouette scores in single-cell data are typically low due to continuous rather than discrete cell states; the metric is used here for relative comparison between resolutions, not as an absolute quality measure.

## Quick Start

```bash
git clone https://github.com/Ekin-Kahraman/single-cell-rnaseq-immune-profiling.git
cd single-cell-rnaseq-immune-profiling
pip install -e .
python run_pipeline.py            # full pipeline (~17s)
python run_pipeline.py --from 4   # resume from step 4
```

## Testing

```bash
pip install -e ".[dev]"
pytest -v
```

7 tests covering QC filtering, normalisation, HVG selection, clustering, and marker gene validation. CI runs on Python 3.10, 3.11, and 3.12.

## Design Decisions

- **Automated annotation** — Clusters are scored against curated PBMC marker gene sets rather than annotated by manual inspection. This makes the pipeline reproducible and removes subjective judgement.
- **Multi-resolution clustering** — Running Leiden at multiple resolutions and picking by silhouette score (with a biological floor) avoids the common problem of choosing an arbitrary resolution.
- **Colourblind-friendly palette** — Okabe-Ito colours throughout.
- **Modular scripts** — Each step is independent. Re-run any step without repeating upstream work.

## Limitations and Future Work

- **No doublet detection.** Scrublet or similar should precede QC in a production pipeline. Omitted here because PBMC 3k is a clean benchmark with negligible doublet rates.
- **No batch correction.** Single-sample dataset. Multi-sample analyses would require Harmony, scVI, or BBKNN.
- **`regress_out` is debatable.** Used here following the original scanpy tutorial, but Luecken & Theis (2019) suggest regression may overcorrect for well-filtered cells. Included for pedagogical alignment with the standard workflow.
- **CD8+ T cells not resolved.** Would require higher clustering resolution or subclustering of the T cell compartment.

## Licence

MIT
