"""Generate an animated 3D UMAP rotation GIF for the README."""

import scanpy as sc
import matplotlib.pyplot as plt
import imageio.v3 as iio
from pathlib import Path
from io import BytesIO

RESULTS_DIR = Path("results")
DOCS_DIR = Path("docs")

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

N_FRAMES = 120
FPS = 24


def compute_3d_umap(adata):
    """Compute 3D UMAP embedding."""
    sc.tl.umap(adata, n_components=3)
    return adata.obsm["X_umap"]


def render_frame(coords, cell_types, azim, elev=25):
    """Render a single frame of the 3D UMAP at a given azimuth angle."""
    fig = plt.figure(figsize=(8, 6), facecolor="white")
    ax = fig.add_subplot(111, projection="3d", facecolor="white")

    for ct in cell_types.cat.categories:
        mask = cell_types == ct
        ax.scatter(
            coords[mask, 0], coords[mask, 1], coords[mask, 2],
            c=PALETTE.get(ct, "#AAAAAA"),
            s=6, alpha=0.8, label=ct, edgecolors="none",
        )

    ax.view_init(elev=elev, azim=azim)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("#EEEEEE")
    ax.yaxis.pane.set_edgecolor("#EEEEEE")
    ax.zaxis.pane.set_edgecolor("#EEEEEE")
    ax.xaxis.line.set_color("#CCCCCC")
    ax.yaxis.line.set_color("#CCCCCC")
    ax.zaxis.line.set_color("#CCCCCC")
    ax.grid(True, alpha=0.15)

    ax.legend(
        loc="upper left", fontsize=7, framealpha=0.7,
        facecolor="white", edgecolor="#DDDDDD",
        markerscale=3,
    )

    ax.set_title("3D UMAP — PBMC Immune Cell Profiling", color="#222222",
                 fontsize=14, fontweight="bold", pad=10)

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    return iio.imread(buf)


def main():
    in_path = RESULTS_DIR / "05_annotated.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    # Need to recompute neighbor graph since preprocessed data is subset
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

    print("Computing 3D UMAP...")
    coords = compute_3d_umap(adata)
    cell_types = adata.obs["cell_type"]

    print(f"Rendering {N_FRAMES} frames...")
    frames = []
    for i in range(N_FRAMES):
        azim = (i / N_FRAMES) * 360
        frame = render_frame(coords, cell_types, azim)
        frames.append(frame)
        if (i + 1) % 30 == 0:
            print(f"  {i + 1}/{N_FRAMES} frames")

    DOCS_DIR.mkdir(exist_ok=True)
    out_path = DOCS_DIR / "umap_3d_rotation.gif"
    print(f"Writing GIF to {out_path}...")
    iio.imwrite(out_path, frames, duration=int(1000 / FPS), loop=0)

    size_mb = out_path.stat().st_size / 1024 / 1024
    print(f"Done! GIF size: {size_mb:.1f} MB")

    return out_path


if __name__ == "__main__":
    main()
