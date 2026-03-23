"""Slow-rotating 3D UMAP for profile banner. Clean, minimal, floaty."""

import scanpy as sc
import matplotlib.pyplot as plt
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

BG = "#0d1117"
N_FRAMES = 120
FPS = 20


def render_frame(coords_3d, cell_types, azim, elev=20):
    fig = plt.figure(figsize=(9, 6), facecolor=BG)
    ax = fig.add_subplot(111, projection="3d", facecolor=BG)

    for ct in cell_types.cat.categories:
        mask = cell_types == ct
        ax.scatter(
            coords_3d[mask, 0], coords_3d[mask, 1], coords_3d[mask, 2],
            c=PALETTE.get(ct, "#AAAAAA"), s=5, alpha=0.8,
            label=ct, edgecolors="none",
        )

    ax.view_init(elev=elev, azim=azim)

    # Clean axes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor(BG)
    ax.yaxis.pane.set_edgecolor(BG)
    ax.zaxis.pane.set_edgecolor(BG)
    ax.xaxis.line.set_color(BG)
    ax.yaxis.line.set_color(BG)
    ax.zaxis.line.set_color(BG)
    ax.grid(False)

    ax.legend(
        loc="upper left", fontsize=7, framealpha=0.0,
        labelcolor="#c9d1d9", edgecolor="none",
        markerscale=3,
    )

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=110, facecolor=BG, edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    return iio.imread(buf)


def main():
    adata = sc.read_h5ad(RESULTS_DIR / "05_annotated.h5ad")
    print(f"Loaded {adata.n_obs} cells")

    # Compute 3D UMAP
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, n_components=3)
    coords_3d = adata.obsm["X_umap"]
    cell_types = adata.obs["cell_type"]

    print(f"Rendering {N_FRAMES} frames...")
    frames = []
    for i in range(N_FRAMES):
        azim = (i / N_FRAMES) * 360
        frame = render_frame(coords_3d, cell_types, azim)
        frames.append(frame)
        if (i + 1) % 30 == 0:
            print(f"  {i + 1}/{N_FRAMES}")

    out_path = Path("docs") / "slow_rotation.gif"
    out_path.parent.mkdir(exist_ok=True)
    iio.imwrite(out_path, frames, duration=int(1000 / FPS), loop=0)
    print(f"Done! {out_path.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
