"""Step 01: Load PBMC 3k dataset and perform quality control."""

import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

# --- Parameters ---
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 3
MAX_PCT_MITO = 5
MAX_GENES_PER_CELL = 2500

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
FIG_DIR = RESULTS_DIR / "figures"


def load_pbmc3k():
    """Download and return the 10X PBMC 3k dataset."""
    DATA_DIR.mkdir(exist_ok=True)
    cache_path = DATA_DIR / "pbmc3k_raw.h5ad"
    if cache_path.exists():
        print(f"Loading cached data from {cache_path}")
        return sc.read_h5ad(cache_path)
    print("Downloading PBMC 3k dataset from 10X Genomics...")
    adata = sc.datasets.pbmc3k()
    adata.write(cache_path)
    print(f"Saved raw data to {cache_path}")
    return adata


def run_qc(adata):
    """Calculate QC metrics and filter cells/genes."""
    print(f"Raw dataset: {adata.n_obs} cells, {adata.n_vars} genes")

    # Ensure gene names are unique
    adata.var_names_make_unique()

    # Annotate mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # Annotate ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )

    # --- QC Plots ---
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    fig.suptitle("Quality Control Metrics (pre-filtering)", fontsize=14, y=1.02)

    axes[0].hist(adata.obs["n_genes_by_counts"], bins=50, color="#2196F3", edgecolor="white")
    axes[0].axvline(MIN_GENES_PER_CELL, color="red", linestyle="--", label=f"min={MIN_GENES_PER_CELL}")
    axes[0].axvline(MAX_GENES_PER_CELL, color="red", linestyle="--", label=f"max={MAX_GENES_PER_CELL}")
    axes[0].set_xlabel("Genes per cell")
    axes[0].set_ylabel("Cells")
    axes[0].legend(fontsize=8)

    axes[1].hist(adata.obs["total_counts"], bins=50, color="#4CAF50", edgecolor="white")
    axes[1].set_xlabel("Total counts per cell")
    axes[1].set_ylabel("Cells")

    axes[2].hist(adata.obs["pct_counts_mt"], bins=50, color="#FF9800", edgecolor="white")
    axes[2].axvline(MAX_PCT_MITO, color="red", linestyle="--", label=f"max={MAX_PCT_MITO}%")
    axes[2].set_xlabel("Mitochondrial %")
    axes[2].set_ylabel("Cells")
    axes[2].legend(fontsize=8)

    axes[3].scatter(
        adata.obs["total_counts"],
        adata.obs["pct_counts_mt"],
        s=1, alpha=0.3, c="#9C27B0", rasterized=True,
    )
    axes[3].axhline(MAX_PCT_MITO, color="red", linestyle="--")
    axes[3].set_xlabel("Total counts")
    axes[3].set_ylabel("Mitochondrial %")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "01_qc_metrics.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved QC plot to {FIG_DIR / '01_qc_metrics.png'}")

    # --- Doublet detection ---
    sc.pp.scrublet(adata, random_state=42)
    n_doublets = adata.obs["predicted_doublet"].sum()
    print(f"Scrublet detected {n_doublets} predicted doublets ({n_doublets/adata.n_obs*100:.1f}%)")

    # --- Filtering ---
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    adata = adata[adata.obs["n_genes_by_counts"] < MAX_GENES_PER_CELL, :].copy()
    adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MITO, :].copy()
    # Remove predicted doublets
    n_doublets_remaining = adata.obs["predicted_doublet"].sum()
    adata = adata[~adata.obs["predicted_doublet"], :].copy()
    n_after = adata.n_obs
    print(f"Removed {n_doublets_remaining} doublets from filtered cells")

    print(f"Filtered: {n_before} -> {n_after} cells ({n_before - n_after} removed)")
    print(f"Genes remaining: {adata.n_vars}")

    return adata


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    adata = load_pbmc3k()
    adata = run_qc(adata)
    out_path = RESULTS_DIR / "01_filtered.h5ad"
    adata.write(out_path)
    print(f"Saved filtered data to {out_path}")
    return adata


if __name__ == "__main__":
    main()
