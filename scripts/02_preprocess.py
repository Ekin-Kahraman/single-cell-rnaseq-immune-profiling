"""Step 02: Normalise, find highly variable genes, and scale."""

import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

# --- Parameters ---
TARGET_SUM = 1e4
N_TOP_GENES = 2000

RESULTS_DIR = Path("results")
FIG_DIR = RESULTS_DIR / "figures"


def preprocess(adata):
    """Normalise, log-transform, select HVGs, regress, and scale."""
    print(f"Input: {adata.n_obs} cells, {adata.n_vars} genes")

    # Normalise to target_sum counts per cell
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    print(f"Normalized to {TARGET_SUM:.0f} counts per cell")

    # Log-transform
    sc.pp.log1p(adata)

    # Store normalized+log data for later use in DE and plotting
    adata.raw = adata

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)

    n_hvg = adata.var["highly_variable"].sum()
    print(f"Selected {n_hvg} highly variable genes")

    # Plot HVG selection
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(FIG_DIR / "02_hvg_selection.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"Saved HVG plot to {FIG_DIR / '02_hvg_selection.png'}")

    # Subset to HVGs
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Regress out unwanted variation
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    print("Regressed out total_counts and pct_counts_mt")

    # Scale to unit variance, clip at max_value=10
    sc.pp.scale(adata, max_value=10)
    print("Scaled data")

    return adata


def main():
    in_path = RESULTS_DIR / "01_filtered.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    adata = preprocess(adata)

    out_path = RESULTS_DIR / "02_preprocessed.h5ad"
    adata.write(out_path)
    print(f"Saved preprocessed data to {out_path}")
    return adata


if __name__ == "__main__":
    main()
