"""Step 06: PAGA trajectory inference and diffusion pseudotime."""

import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS_DIR = Path("results")
FIG_DIR = RESULTS_DIR / "figures"


def compute_trajectory(adata):
    """Run PAGA and diffusion pseudotime on annotated data."""
    # PAGA — partition-based graph abstraction
    sc.tl.paga(adata, groups="cell_type")
    sc.pl.paga(adata, threshold=0.03, show=False)
    plt.savefig(FIG_DIR / "06_paga_graph.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIG_DIR / "06_paga_graph.pdf", bbox_inches="tight")
    plt.close()
    print("Computed PAGA graph")

    # Reinitialise UMAP using PAGA as initialisation for cleaner layout
    sc.tl.umap(adata, init_pos="paga", random_state=42)

    # Diffusion pseudotime — root in CD14+ monocytes (most primitive in PBMC)
    sc.tl.diffmap(adata)

    # Find a root cell in CD14+ monocytes
    monocyte_mask = adata.obs["cell_type"] == "CD14+ Monocytes"
    if monocyte_mask.any():
        adata.uns["iroot"] = monocyte_mask.values.nonzero()[0][0]
        sc.tl.dpt(adata)
        print("Computed diffusion pseudotime (rooted in CD14+ monocytes)")
    else:
        print("Warning: No CD14+ monocytes found, skipping diffusion pseudotime")

    return adata


def plot_trajectory(adata):
    """Plot PAGA-initialised UMAP and pseudotime."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # PAGA-initialised UMAP coloured by cell type
    sc.pl.umap(adata, color="cell_type", legend_loc="right margin",
               frameon=False, ax=axes[0], show=False, title="Cell types (PAGA layout)")

    # PAGA connectivity overlaid on UMAP
    sc.pl.paga(adata, pos=adata.uns["paga"]["pos"], threshold=0.03,
               node_size_scale=1.5, ax=axes[1], show=False,
               title="PAGA connectivity")

    # Diffusion pseudotime
    if "dpt_pseudotime" in adata.obs:
        sc.pl.umap(adata, color="dpt_pseudotime", frameon=False,
                   ax=axes[2], show=False, title="Diffusion pseudotime",
                   color_map="viridis")
    else:
        axes[2].text(0.5, 0.5, "Pseudotime not computed",
                     ha="center", va="center", transform=axes[2].transAxes)
        axes[2].set_axis_off()

    fig.tight_layout()
    fig.savefig(FIG_DIR / "06_trajectory.png", dpi=150, bbox_inches="tight")
    fig.savefig(FIG_DIR / "06_trajectory.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved trajectory plots to {FIG_DIR / '06_trajectory.png'}")


def main():
    in_path = RESULTS_DIR / "05_annotated.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    adata = compute_trajectory(adata)
    plot_trajectory(adata)

    out_path = RESULTS_DIR / "06_trajectory.h5ad"
    adata.write(out_path)
    print(f"Saved trajectory data to {out_path}")
    return adata


if __name__ == "__main__":
    main()
