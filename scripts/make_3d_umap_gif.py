"""Generate an animated 3D UMAP rotation GIF for the README."""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import imageio.v3 as iio
from pathlib import Path
from io import BytesIO
from palette import PALETTE

RESULTS_DIR = Path("results")
DOCS_DIR = Path("docs")

N_FRAMES = 120
FPS = 24


def compute_3d_umap(adata):
    """Compute 3D UMAP embedding with fixed seed."""
    sc.tl.umap(adata, n_components=3, random_state=42)
    return adata.obsm["X_umap"]


def render_frame(coords, cell_types, azim, elev):
    """Render a single frame of the 3D UMAP."""
    fig = plt.figure(figsize=(7, 7), facecolor="#0D1117")
    ax = fig.add_subplot(111, projection="3d", facecolor="#0D1117")

    # Sort cell types so smaller populations render on top
    type_order = cell_types.value_counts().index[::-1]

    for ct in type_order:
        mask = cell_types == ct
        colour = PALETTE.get(ct, "#AAAAAA")
        ax.scatter(
            coords[mask, 0], coords[mask, 1], coords[mask, 2],
            c=colour, s=10, alpha=0.85, label=ct,
            edgecolors="none", rasterized=True, depthshade=True,
        )

    ax.view_init(elev=elev, azim=azim)

    # Clean axes — no ticks, no panes, no grid
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("#0D1117")
    ax.yaxis.pane.set_edgecolor("#0D1117")
    ax.zaxis.pane.set_edgecolor("#0D1117")
    ax.xaxis.line.set_color("#0D1117")
    ax.yaxis.line.set_color("#0D1117")
    ax.zaxis.line.set_color("#0D1117")
    ax.grid(False)

    # Legend — bottom right, out of the way
    legend = ax.legend(
        loc="lower right", fontsize=8, framealpha=0.85,
        facecolor="#161B22", edgecolor="#30363D",
        labelcolor="white", markerscale=3,
    )
    legend.get_frame().set_linewidth(0.5)

    ax.set_title("PBMC Immune Cell Profiling — 3D UMAP",
                 color="white", fontsize=13, fontweight="bold", pad=15)

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight",
                facecolor="#0D1117", edgecolor="none", pad_inches=0.3)
    plt.close(fig)
    buf.seek(0)
    return iio.imread(buf)


def main():
    in_path = RESULTS_DIR / "05_annotated.h5ad"
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {in_path}")

    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

    print("Computing 3D UMAP...")
    coords = compute_3d_umap(adata)
    cell_types = adata.obs["cell_type"]

    print(f"Rendering {N_FRAMES} frames...")
    frames = []
    for i in range(N_FRAMES):
        azim = (i / N_FRAMES) * 360
        # Gentle elevation oscillation: 20° to 35° and back
        elev = 27.5 + 7.5 * np.sin(2 * np.pi * i / N_FRAMES)
        frame = render_frame(coords, cell_types, azim, elev)
        frames.append(frame)
        if (i + 1) % 45 == 0:
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
