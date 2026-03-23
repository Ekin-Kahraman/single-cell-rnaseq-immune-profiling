"""Generate an animated profile banner: UMAP atlas reveal + volcano plot."""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import imageio.v3 as iio
from pathlib import Path
from io import BytesIO

RESULTS_DIR = Path("results")

PALETTE = {
    "CD4+ T cells": "#E69F00",
    "CD8+ T cells": "#56B4E9",
    "NK cells": "#009E73",
    "B cells": "#F0E442",
    "CD14+ Monocytes": "#0072B2",
    "FCGR3A+ Monocytes": "#D55E00",
    "Dendritic cells": "#CC79A7",
    "Megakaryocytes": "#999999",
}

BG = "#0d1117"  # GitHub dark mode background
TEXT = "#e6edf3"
GRID = "#21262d"

FPS = 20


def render_frame(umap_coords, cell_types, categories, reveal_progress,
                 xlim=None, ylim=None):
    """Render one frame of the atlas reveal."""
    fig, ax = plt.subplots(figsize=(10, 7), facecolor=BG)
    ax.set_facecolor(BG)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    n_types = len(categories)
    # Each cluster gets an equal slice of the progress bar
    cluster_slice = 1.0 / n_types

    for i, ct in enumerate(categories):
        mask = cell_types == ct
        coords = umap_coords[mask]
        colour = PALETTE.get(ct, "#AAAAAA")

        # Calculate this cluster's visibility (0 to 1)
        cluster_start = i * cluster_slice
        if reveal_progress < cluster_start:
            continue  # not yet visible

        # Fade in: 0 at cluster_start, 1 at cluster_end
        alpha = min(1.0, (reveal_progress - cluster_start) / cluster_slice)
        alpha = alpha ** 0.5  # ease-in curve

        ax.scatter(
            coords[:, 0], coords[:, 1],
            c=colour, s=5, alpha=alpha * 0.85,
            edgecolors="none", rasterized=True,
        )

        # Show label once cluster is >60% visible
        if alpha > 0.6:
            cx, cy = coords[:, 0].mean(), coords[:, 1].mean()
            ax.text(
                cx, cy, ct, fontsize=8, color="white",
                ha="center", va="center", fontweight="bold",
                alpha=min(1.0, (alpha - 0.6) / 0.4),
                path_effects=[
                    pe.withStroke(linewidth=2.5, foreground=BG, alpha=min(1.0, (alpha - 0.6) / 0.4))
                ],
            )

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Title with fade
    title_alpha = min(1.0, reveal_progress * 3)  # fades in early
    ax.set_title(
        "PBMC Immune Cell Atlas  ·  2,638 cells  ·  6 cell types",
        color=TEXT, fontsize=13, fontweight="bold", pad=15,
        alpha=title_alpha,
    )

    # Subtitle
    if reveal_progress > 0.1:
        sub_alpha = min(1.0, (reveal_progress - 0.1) * 5)
        ax.text(
            0.5, -0.02, "scanpy · Leiden clustering · automated marker-based annotation",
            transform=ax.transAxes, ha="center", fontsize=9,
            color=TEXT, alpha=sub_alpha * 0.6,
        )

    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.05)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120, facecolor=BG,
                edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    return iio.imread(buf)


def main():
    adata = sc.read_h5ad(RESULTS_DIR / "05_annotated.h5ad")
    print(f"Loaded {adata.n_obs} cells")

    umap_coords = adata.obsm["X_umap"]
    cell_types = adata.obs["cell_type"]
    categories = cell_types.cat.categories.tolist()

    # Sort categories by size (largest first) for dramatic reveal
    sizes = cell_types.value_counts()
    categories = sizes.index.tolist()

    # Precompute axis limits so they're consistent across frames
    pad = 1.5
    xlim = (umap_coords[:, 0].min() - pad, umap_coords[:, 0].max() + pad)
    ylim = (umap_coords[:, 1].min() - pad, umap_coords[:, 1].max() + pad)

    # Build frames: reveal phase + hold phase
    n_reveal = 80
    n_hold = 40
    total = n_reveal + n_hold

    print(f"Rendering {total} frames...")
    frames = []
    for i in range(total):
        if i < n_reveal:
            # Start at 0.05 so first frame already has visible cells
            progress = 0.05 + 0.95 * (i / (n_reveal - 1))
        else:
            progress = 1.0

        frame = render_frame(umap_coords, cell_types, categories, progress,
                             xlim=xlim, ylim=ylim)
        frames.append(frame)
        if (i + 1) % 30 == 0:
            print(f"  {i + 1}/{total}")

    out_path = Path("docs") / "profile_banner.gif"
    out_path.parent.mkdir(exist_ok=True)
    print(f"Writing GIF to {out_path}...")
    iio.imwrite(out_path, frames, duration=int(1000 / FPS), loop=0)

    size_mb = out_path.stat().st_size / 1024 / 1024
    print(f"Done! {size_mb:.1f} MB")


if __name__ == "__main__":
    main()
