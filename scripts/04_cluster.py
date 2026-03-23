"""Step 04: Leiden clustering at multiple resolutions with evaluation."""

import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from pathlib import Path

# --- Parameters ---
RESOLUTIONS = [0.3, 0.5, 0.8, 1.0, 1.2]

RESULTS_DIR = Path("results")
FIG_DIR = RESULTS_DIR / "figures"


def cluster_multi_resolution(adata):
    """Run Leiden clustering at multiple resolutions and pick the best."""
    scores = {}

    for res in RESOLUTIONS:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        n_clusters = adata.obs[key].nunique()

        if n_clusters > 1:
            sil = silhouette_score(
                adata.obsm["X_pca"][:, : adata.uns.get("n_pcs_selected", 40)],
                adata.obs[key],
                sample_size=min(5000, adata.n_obs),
                random_state=42,
            )
        else:
            sil = -1.0

        scores[res] = {"silhouette": sil, "n_clusters": n_clusters}
        print(f"  Resolution {res}: {n_clusters} clusters, silhouette={sil:.3f}")

    # Select best resolution: prefer silhouette but require >= 5 clusters
    # (PBMCs have at least 5 distinct cell types)
    candidates = {r: s for r, s in scores.items() if s["n_clusters"] >= 5}
    if not candidates:
        candidates = scores
    best_res = max(candidates, key=lambda r: candidates[r]["silhouette"])
    best_key = f"leiden_{best_res}"
    adata.obs["leiden"] = adata.obs[best_key].copy()
    adata.uns["clustering"] = {
        "best_resolution": best_res,
        "scores": {str(k): v for k, v in scores.items()},
    }
    print(f"\nBest resolution: {best_res} ({scores[best_res]['n_clusters']} clusters, "
          f"silhouette={scores[best_res]['silhouette']:.3f})")

    return adata, scores


def plot_clustering(adata, scores):
    """Plot clustering results: resolution comparison + UMAP."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Silhouette scores by resolution
    resolutions = list(scores.keys())
    sil_scores = [scores[r]["silhouette"] for r in resolutions]
    n_clusters = [scores[r]["n_clusters"] for r in resolutions]
    best_res = adata.uns["clustering"]["best_resolution"]

    colors = ["#E53935" if r == best_res else "#90CAF9" for r in resolutions]
    axes[0].bar([str(r) for r in resolutions], sil_scores, color=colors, edgecolor="white")
    axes[0].set_xlabel("Resolution")
    axes[0].set_ylabel("Silhouette Score")
    axes[0].set_title("Clustering Quality by Resolution")

    # Number of clusters by resolution
    axes[1].plot(resolutions, n_clusters, "o-", color="#4CAF50", linewidth=2, markersize=8)
    axes[1].axvline(best_res, color="red", linestyle="--", alpha=0.5)
    axes[1].set_xlabel("Resolution")
    axes[1].set_ylabel("Number of Clusters")
    axes[1].set_title("Cluster Count by Resolution")

    # UMAP with best clustering
    sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontsize=10,
               legend_fontoutline=2, frameon=False, ax=axes[2], show=False,
               title=f"Leiden (res={best_res})")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "04_clustering.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved clustering plot to {FIG_DIR / '04_clustering.png'}")


def main():
    in_path = RESULTS_DIR / "03_reduced.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    adata, scores = cluster_multi_resolution(adata)
    plot_clustering(adata, scores)

    out_path = RESULTS_DIR / "04_clustered.h5ad"
    adata.write(out_path)
    print(f"Saved clustered data to {out_path}")
    return adata


if __name__ == "__main__":
    main()
