"""Display labels that distinguish broad annotation from T-cell refinement."""

import pandas as pd

T_CELL_LABELS = {"CD4+ T cells", "CD8+ T cells"}
SUBTYPE_LABELS = {"CD4+ T": "CD4-like T cells", "CD8+ T": "CD8-like T cells"}


def broad_cell_types(adata):
    labels = adata.obs["cell_type"].astype(str).copy()
    return labels.mask(labels.isin(T_CELL_LABELS), "T cells").astype("category")


def refined_cell_types(adata, t_cells):
    """Join subtype evidence by cell ID, rejecting missing or unrelated cells."""
    expected = adata.obs_names[adata.obs["cell_type"].isin(T_CELL_LABELS)]
    if (not adata.obs_names.is_unique or not t_cells.obs_names.is_unique
            or set(expected) != set(t_cells.obs_names)):
        raise ValueError("T-cell refinement must match every broad T cell exactly once.")
    if "t_subtype" not in t_cells.obs:
        raise ValueError("T-cell refinement is missing t_subtype labels.")
    subtypes = t_cells.obs["t_subtype"].astype(str).map(SUBTYPE_LABELS)
    if subtypes.isna().any():
        raise ValueError("T-cell refinement contains missing or unsupported subtype labels.")
    labels = broad_cell_types(adata).astype(str)
    labels.loc[subtypes.index] = subtypes
    return pd.Series(labels, index=adata.obs_names, name="display_cell_type").astype("category")
