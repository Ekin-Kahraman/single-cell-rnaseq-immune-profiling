# Single-Cell RNA-seq Immune Cell Profiling

End-to-end single-cell RNA-seq analysis pipeline in Python using [scanpy](https://scanpy.readthedocs.io/). Demonstrates quality control, normalization, dimensionality reduction, clustering with automated resolution selection, and marker-based cell type annotation on human PBMC data.

![Publication Figure](results/figures/06_publication_figure.png)

## Dataset

**10X Genomics PBMC 3k** -- 2,700 peripheral blood mononuclear cells from a healthy donor, sequenced on the Chromium platform. A standard benchmark dataset for single-cell analysis pipelines.

- **Source**: [10X Genomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)
- **Reference**: Zheng et al. (2017) *Nature Communications*

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 01 | `scripts/01_load_and_qc.py` | Download data, calculate QC metrics (genes/cell, counts, mito %), filter |
| 02 | `scripts/02_preprocess.py` | Normalize (10k/cell), log-transform, select 2,000 HVGs, regress covariates, scale |
| 03 | `scripts/03_reduce_dimensions.py` | PCA (40 components), neighbor graph, UMAP embedding |
| 04 | `scripts/04_cluster.py` | Leiden clustering at 5 resolutions, silhouette-based selection (min 5 clusters) |
| 05 | `scripts/05_annotate_cell_types.py` | Wilcoxon DE, marker gene scoring, automated cell type assignment |
| 06 | `scripts/06_publication_figures.py` | Multi-panel publication figure (UMAP, composition, heatmap) |

## Results

| Cell Type | Count | Proportion |
|-----------|-------|------------|
| CD4+ T cells | 1,195 | 45.3% |
| CD14+ Monocytes | 464 | 17.6% |
| NK cells | 419 | 15.9% |
| B cells | 342 | 13.0% |
| FCGR3A+ Monocytes | 180 | 6.8% |
| Dendritic cells | 38 | 1.4% |

Best clustering: Leiden resolution 0.5, 6 clusters, silhouette score 0.196.

## Quick Start

```bash
# Clone
git clone https://github.com/Ekin-Kahraman/single-cell-rnaseq-immune-profiling.git
cd single-cell-rnaseq-immune-profiling

# Install
pip install -e .

# Run full pipeline (~17 seconds)
python run_pipeline.py

# Resume from a specific step
python run_pipeline.py --from 4
```

## Requirements

- Python >= 3.10
- scanpy >= 1.10
- leidenalg >= 0.10
- scikit-learn >= 1.3

Full dependency list in `pyproject.toml`.

## Testing

```bash
pip install -e ".[dev]"
pytest -v
```

## Project Structure

```
.
├── scripts/           # Analysis pipeline (01-06)
├── tests/             # pytest test suite
├── data/              # Raw data (auto-downloaded)
├── results/           # Output: h5ad files, figures, CSVs
├── run_pipeline.py    # Pipeline orchestrator
├── pyproject.toml     # Dependencies and project config
└── CITATION.cff       # Citation metadata
```

## Key Design Decisions

- **Automated cell type annotation**: Clusters are assigned to cell types by scoring against curated PBMC marker gene sets, not manual inspection.
- **Multi-resolution clustering**: Leiden is run at 5 resolutions (0.3-1.2) and the best is selected by silhouette score with a biological floor of 5 clusters.
- **Colorblind-friendly palette**: Publication figures use the Okabe-Ito palette for accessibility.
- **Modular scripts**: Each step reads the previous step's output from disk. Steps can be re-run independently.

## License

MIT
