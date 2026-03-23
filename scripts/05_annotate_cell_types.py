"""Step 05: Marker gene identification and automated cell type annotation."""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# --- Parameters ---
N_MARKER_GENES = 25
MIN_FOLD_CHANGE = 1.5

RESULTS_DIR = Path("results")
FIG_DIR = RESULTS_DIR / "figures"

# Known PBMC marker gene sets for automated annotation
PBMC_MARKERS = {
    "CD4+ T cells": ["CD3D", "CD3E", "IL7R", "CD4"],
    "CD8+ T cells": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMK"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "KLRB1"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "CD14+ Monocytes": ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7", "LST1"],
    "Dendritic cells": ["FCER1A", "CST3", "IL3RA"],
    "Megakaryocytes": ["PPBP", "PF4", "GP9"],
}


def find_markers(adata):
    """Identify marker genes for each cluster using Wilcoxon test."""
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True)

    # Extract top markers per cluster
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers = markers[markers["logfoldchanges"] > np.log2(MIN_FOLD_CHANGE)]
    markers = markers[markers["pvals_adj"] < 0.05]

    n_per_cluster = markers.groupby("group").size()
    print("Significant markers per cluster (log2FC > {:.1f}, padj < 0.05):".format(
        np.log2(MIN_FOLD_CHANGE)))
    for cluster, count in n_per_cluster.items():
        print(f"  Cluster {cluster}: {count} genes")

    # Save full marker table
    markers.to_csv(RESULTS_DIR / "marker_genes.csv", index=False)
    print(f"Saved marker gene table to {RESULTS_DIR / 'marker_genes.csv'}")

    return markers


def score_cell_types(adata):
    """Score each cell for known cell type signatures and assign labels."""
    available_genes = set(adata.raw.var_names) if adata.raw else set(adata.var_names)

    scores = {}
    for cell_type, genes in PBMC_MARKERS.items():
        present = [g for g in genes if g in available_genes]
        if len(present) < 2:
            print(f"  Skipping {cell_type}: only {len(present)} markers found")
            continue
        score_name = cell_type.replace("+", "p").replace(" ", "_").lower()
        sc.tl.score_genes(adata, gene_list=present, score_name=score_name, use_raw=True)
        scores[cell_type] = score_name
        print(f"  Scored {cell_type} ({len(present)}/{len(genes)} markers)")

    # Assign each cluster to the cell type with highest mean score
    cluster_labels = {}
    for cluster in adata.obs["leiden"].cat.categories:
        mask = adata.obs["leiden"] == cluster
        best_type = None
        best_score = -np.inf
        for cell_type, score_name in scores.items():
            mean_score = adata.obs.loc[mask, score_name].mean()
            if mean_score > best_score:
                best_score = mean_score
                best_type = cell_type
        cluster_labels[cluster] = best_type
        print(f"  Cluster {cluster} -> {best_type} (score={best_score:.3f})")

    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_labels).astype("category")

    # Summary
    print("\nCell type composition:")
    composition = adata.obs["cell_type"].value_counts()
    for ct, count in composition.items():
        pct = count / adata.n_obs * 100
        print(f"  {ct}: {count} cells ({pct:.1f}%)")

    composition.to_csv(RESULTS_DIR / "cell_type_composition.csv")
    return adata


def plot_annotation(adata):
    """Plot annotation results: markers dotplot + UMAP."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Flatten marker genes for dotplot (keep only those present)
    available = set(adata.raw.var_names) if adata.raw else set(adata.var_names)
    marker_dict = {}
    for ct, genes in PBMC_MARKERS.items():
        present = [g for g in genes if g in available]
        if present:
            marker_dict[ct] = present

    # Dotplot of canonical markers by cluster
    fig_dot = sc.pl.dotplot(
        adata, var_names=marker_dict, groupby="leiden",
        standard_scale="var", use_raw=True, return_fig=True,
    )
    fig_dot.savefig(FIG_DIR / "05_marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"Saved marker dotplot to {FIG_DIR / '05_marker_dotplot.png'}")

    # UMAP colored by cell type
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2,
               frameon=False, ax=axes[0], show=False, title="Clusters")
    sc.pl.umap(adata, color="cell_type", legend_loc="right margin",
               frameon=False, ax=axes[1], show=False, title="Cell Types")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "05_cell_types_umap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved cell type UMAP to {FIG_DIR / '05_cell_types_umap.png'}")


def main():
    in_path = RESULTS_DIR / "04_clustered.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    print("\n--- Finding marker genes ---")
    find_markers(adata)

    print("\n--- Scoring and annotating cell types ---")
    adata = score_cell_types(adata)

    print("\n--- Generating annotation plots ---")
    plot_annotation(adata)

    out_path = RESULTS_DIR / "05_annotated.h5ad"
    adata.write(out_path)
    print(f"\nSaved annotated data to {out_path}")
    return adata


if __name__ == "__main__":
    main()
