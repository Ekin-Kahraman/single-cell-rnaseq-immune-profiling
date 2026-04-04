"""Step 03: PCA, neighbor graph, and UMAP embedding."""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# --- Parameters ---
N_PCS = 40
N_NEIGHBORS = 15

RESULTS_DIR = Path("results")
FIG_DIR = RESULTS_DIR / "figures"


def reduce_dimensions(adata):
    """Compute PCA, build neighbor graph, and embed with UMAP."""
    n_comps = min(N_PCS, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
    print(f"Computed PCA with {n_comps} components")

    # Plot variance ratio
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    variance_ratio = adata.uns["pca"]["variance_ratio"]
    cumulative = np.cumsum(variance_ratio)

    axes[0].bar(range(1, len(variance_ratio) + 1), variance_ratio, color="#2196F3", edgecolor="white")
    axes[0].axvline(N_PCS, color="red", linestyle="--", label=f"n_pcs={N_PCS}")
    axes[0].set_xlabel("Principal Component")
    axes[0].set_ylabel("Variance Ratio")
    axes[0].set_title("Scree Plot")
    axes[0].legend()

    axes[1].plot(range(1, len(cumulative) + 1), cumulative, color="#4CAF50", linewidth=2)
    axes[1].axhline(0.90, color="gray", linestyle=":", label="90% threshold")
    axes[1].axvline(N_PCS, color="red", linestyle="--", label=f"n_pcs={N_PCS}")
    axes[1].set_xlabel("Principal Component")
    axes[1].set_ylabel("Cumulative Variance")
    axes[1].set_title("Cumulative Variance Explained")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_pca_variance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved PCA variance plot to {FIG_DIR / '03_pca_variance.png'}")

    # Build neighbor graph
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    print(f"Built neighbor graph (n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS})")

    # UMAP embedding (seed for reproducibility)
    sc.tl.umap(adata, random_state=42)
    print("Computed UMAP embedding")

    # Store n_pcs used for downstream reference
    adata.uns["n_pcs_selected"] = N_PCS

    return adata


def main():
    in_path = RESULTS_DIR / "02_preprocessed.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    adata = reduce_dimensions(adata)

    out_path = RESULTS_DIR / "03_reduced.h5ad"
    adata.write(out_path)
    print(f"Saved reduced data to {out_path}")
    return adata


if __name__ == "__main__":
    main()
