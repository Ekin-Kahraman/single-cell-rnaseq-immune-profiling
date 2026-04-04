"""Step 07: Subcluster the T cell compartment to resolve CD4+/CD8+ populations."""

import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS_DIR = Path("results")
FIG_DIR = RESULTS_DIR / "figures"

# Markers for T cell subtypes
T_CELL_MARKERS = {
    "CD4+ T": ["IL7R", "CD4", "TCF7", "LEF1"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "NKG7"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4"],
    "Naive T": ["CCR7", "SELL", "TCF7"],
    "Memory T": ["S100A4", "ANXA1", "IL7R"],
}


def subcluster_t_cells(adata):
    """Extract T cell clusters and subcluster to resolve CD4+/CD8+."""
    # Extract T cells (currently labelled as CD4+ T cells — CD8+ merged in)
    t_mask = adata.obs["cell_type"] == "CD4+ T cells"
    adata_t = adata[t_mask].copy()
    print(f"Extracted {adata_t.n_obs} T cells for subclustering")

    # Recompute neighbours and UMAP on T cell subset
    sc.pp.neighbors(adata_t, n_neighbors=15, n_pcs=20)
    sc.tl.umap(adata_t, random_state=42)
    sc.tl.leiden(adata_t, resolution=0.5, key_added="t_subcluster",
                 flavor="igraph", n_iterations=2, directed=False, random_state=42)

    # Score for CD4 vs CD8
    available = set(adata_t.raw.var_names) if adata_t.raw else set(adata_t.var_names)
    cd4_genes = [g for g in ["IL7R", "CD4", "TCF7", "LEF1"] if g in available]
    cd8_genes = [g for g in ["CD8A", "CD8B", "GZMK", "GZMA"] if g in available]

    sc.tl.score_genes(adata_t, gene_list=cd4_genes, score_name="cd4_score", use_raw=True)
    sc.tl.score_genes(adata_t, gene_list=cd8_genes, score_name="cd8_score", use_raw=True)

    # Assign subtype
    adata_t.obs["t_subtype"] = "CD4+ T"
    cd8_mask = adata_t.obs["cd8_score"] > adata_t.obs["cd4_score"]
    adata_t.obs.loc[cd8_mask, "t_subtype"] = "CD8+ T"

    n_cd4 = (adata_t.obs["t_subtype"] == "CD4+ T").sum()
    n_cd8 = (adata_t.obs["t_subtype"] == "CD8+ T").sum()
    print(f"  CD4+ T: {n_cd4} cells ({n_cd4/adata_t.n_obs*100:.1f}%)")
    print(f"  CD8+ T: {n_cd8} cells ({n_cd8/adata_t.n_obs*100:.1f}%)")

    return adata_t


def plot_subclustering(adata_t):
    """Plot T cell subclustering results."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    sc.pl.umap(adata_t, color="t_subtype", frameon=False,
               ax=axes[0], show=False, title="T cell subtypes")

    sc.pl.umap(adata_t, color="cd4_score", frameon=False,
               ax=axes[1], show=False, title="CD4 score", color_map="Reds")

    sc.pl.umap(adata_t, color="cd8_score", frameon=False,
               ax=axes[2], show=False, title="CD8 score", color_map="Blues")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "07_t_cell_subclustering.png", dpi=150, bbox_inches="tight")
    fig.savefig(FIG_DIR / "07_t_cell_subclustering.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved T cell subclustering plots to {FIG_DIR}")

    # Marker dotplot for T subtypes
    available = set(adata_t.raw.var_names) if adata_t.raw else set(adata_t.var_names)
    markers_present = {}
    for subtype, genes in T_CELL_MARKERS.items():
        present = [g for g in genes if g in available]
        if present:
            markers_present[subtype] = present

    if markers_present:
        fig_dot = sc.pl.dotplot(adata_t, var_names=markers_present,
                                groupby="t_subtype", standard_scale="var",
                                use_raw=True, return_fig=True)
        fig_dot.savefig(FIG_DIR / "07_t_cell_markers.png", dpi=150, bbox_inches="tight")
        plt.close("all")
        print("Saved T cell marker dotplot")


def main():
    in_path = RESULTS_DIR / "06_trajectory.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    adata_t = subcluster_t_cells(adata)
    plot_subclustering(adata_t)

    out_path = RESULTS_DIR / "07_t_cells.h5ad"
    adata_t.write(out_path)
    print(f"Saved T cell data to {out_path}")
    return adata_t


if __name__ == "__main__":
    main()
