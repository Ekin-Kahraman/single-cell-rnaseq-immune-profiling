"""Step 08: Generate publication-quality multi-panel figures."""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from pathlib import Path

from palette import PALETTE

RESULTS_DIR = Path("results")
FIG_DIR = RESULTS_DIR / "figures"

# Key marker genes for heatmap
KEY_MARKERS = [
    "CD3D", "IL7R", "CD8A", "NKG7", "GNLY",
    "MS4A1", "CD79A", "CD14", "LYZ", "FCGR3A",
    "FCER1A", "CST3", "PPBP",
]


def make_figure(adata):
    """Create a multi-panel publication figure."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.45)

    # A: UMAP by cluster
    ax_a = fig.add_subplot(gs[0, 0])
    sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontsize=8,
               legend_fontoutline=2, frameon=False, ax=ax_a, show=False, title="")
    ax_a.set_title("A  Leiden Clusters", fontsize=12, fontweight="bold", loc="left")

    # B: UMAP by cell type (legend on data to avoid overlap with panel C)
    ax_b = fig.add_subplot(gs[0, 1])
    present_types = adata.obs["cell_type"].cat.categories.tolist()
    palette = [PALETTE.get(ct, "#AAAAAA") for ct in present_types]
    sc.pl.umap(adata, color="cell_type", palette=palette,
               legend_loc="on data", legend_fontsize=6, legend_fontoutline=2,
               frameon=False, ax=ax_b, show=False, title="")
    ax_b.set_title("B  Cell Types", fontsize=12, fontweight="bold", loc="left")

    # C: Cell type composition bar chart
    ax_c = fig.add_subplot(gs[0, 2])
    composition = adata.obs["cell_type"].value_counts()
    colors = [PALETTE.get(ct, "#AAAAAA") for ct in composition.index]
    bars = ax_c.barh(range(len(composition)), composition.values, color=colors, edgecolor="white")
    ax_c.set_yticks(range(len(composition)))
    ax_c.set_yticklabels(composition.index, fontsize=8)
    ax_c.set_xlabel("Number of Cells")
    ax_c.set_title("C  Cell Type Composition", fontsize=12, fontweight="bold", loc="left")
    ax_c.invert_yaxis()
    for bar, count in zip(bars, composition.values):
        ax_c.text(bar.get_width() + 5, bar.get_y() + bar.get_height() / 2,
                  str(count), va="center", fontsize=8)

    # D: Marker gene heatmap (bottom row, spanning 2 columns)
    ax_d = fig.add_subplot(gs[1, :2])
    available = set(adata.raw.var_names) if adata.raw else set(adata.var_names)
    markers_present = [g for g in KEY_MARKERS if g in available]

    if markers_present and adata.raw is not None:
        # Build expression matrix: mean expression per cell type
        raw_df = pd.DataFrame(
            adata.raw[:, markers_present].X.toarray()
            if hasattr(adata.raw[:, markers_present].X, "toarray")
            else adata.raw[:, markers_present].X,
            index=adata.obs.index,
            columns=markers_present,
        )
        raw_df["cell_type"] = adata.obs["cell_type"].values
        mean_expr = raw_df.groupby("cell_type")[markers_present].mean()

        # Z-score per gene for visualisation
        from scipy.stats import zscore
        z_expr = mean_expr.apply(zscore, axis=0)

        im = ax_d.imshow(z_expr.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax_d.set_xticks(range(len(markers_present)))
        ax_d.set_xticklabels(markers_present, rotation=45, ha="right", fontsize=9)
        ax_d.set_yticks(range(len(mean_expr)))
        ax_d.set_yticklabels(mean_expr.index, fontsize=8)
        plt.colorbar(im, ax=ax_d, label="Z-score", shrink=0.8)
    ax_d.set_title("D  Marker Gene Expression", fontsize=12, fontweight="bold", loc="left")

    # E: QC summary
    ax_e = fig.add_subplot(gs[1, 2])
    ax_e.axis("off")
    stats = [
        ("Cells", f"{adata.n_obs:,}"),
        ("Genes (HVG)", f"{adata.n_vars:,}"),
        ("Clusters", str(adata.obs["leiden"].nunique())),
        ("Cell types", str(adata.obs["cell_type"].nunique())),
    ]
    if "clustering" in adata.uns:
        stats.append(("Resolution", str(adata.uns["clustering"]["best_resolution"])))
        sil = adata.uns["clustering"]["scores"][
            str(adata.uns["clustering"]["best_resolution"])
        ]["silhouette"]
        stats.append(("Silhouette", f"{sil:.3f}"))

    ax_e.set_title("E  Pipeline Summary", fontsize=12, fontweight="bold", loc="left")
    for i, (label, value) in enumerate(stats):
        y = 0.85 - i * 0.12
        ax_e.text(0.1, y, label, fontsize=11, fontweight="bold", transform=ax_e.transAxes)
        ax_e.text(0.6, y, value, fontsize=11, transform=ax_e.transAxes)

    fig.savefig(FIG_DIR / "08_publication_figure.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "08_publication_figure.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved publication figure to {FIG_DIR / '08_publication_figure.png'}")
    print(f"Saved publication figure to {FIG_DIR / '08_publication_figure.pdf'}")


def main():
    in_path = RESULTS_DIR / "06_trajectory.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    make_figure(adata)
    return adata


if __name__ == "__main__":
    main()
