#!/usr/bin/env python3
"""
reconstruct_h5ad.py

Assemble an AnnData (.h5ad) from the on-disk dump produced by dump_seurat.R.
Every matrix is memory-mapped, so they are never all resident at once, and the
CSC (genes x cells) -> AnnData (cells x genes) transpose is zero-copy.

Usage: python reconstruct_h5ad.py <dump_dir> <output.h5ad>
"""

import os
import sys
import json
import numpy as np
import scipy.sparse as sp
import anndata as ad
import pandas as pd


def _p(dump_dir, tag, ext):
    return os.path.join(dump_dir, tag + ext)


def load_matrix(meta, dump_dir):
    """Return a scipy sparse (CSC) or numpy array, backed by memmaps."""
    fmt = meta["format"]
    shape = (int(meta["shape"][0]), int(meta["shape"][1]))
    tag = meta["tag"]
    if fmt == "csc":
        x = np.memmap(_p(dump_dir, tag, ".x.f64"), dtype="<f8", mode="r")
        i = np.memmap(_p(dump_dir, tag, ".i.i32"), dtype="<i4", mode="r")
        p = np.memmap(_p(dump_dir, tag, ".p.i32"), dtype="<i4", mode="r")
        return sp.csc_matrix((x, i.astype(np.int32, copy=False),
                              p.astype(np.int32, copy=False)), shape=shape)
    elif fmt == "dense":
        d = np.memmap(_p(dump_dir, tag, ".d.f64"), dtype="<f8", mode="r")
        order = meta.get("order", "F")
        return np.asarray(d).reshape(shape, order=order)
    raise ValueError(f"Unknown matrix format: {fmt}")


def to_cells_by_genes(m):
    """Transpose a (genes x cells) matrix to (cells x genes), zero-copy."""
    if sp.issparse(m):
        # CSC (genes x cells).T -> CSR (cells x genes), shares buffers.
        return m.T.tocsr()
    return np.ascontiguousarray(m.T)


def build_dataframe(spec, index, n):
    def as_list(x):
        # jsonlite auto_unbox collapses length-1 arrays to scalars; re-wrap.
        if x is None:
            return []
        if isinstance(x, (list, tuple)):
            return list(x)
        return [x]

    # jsonlite serializes an empty named list as [] (a list), not {}.
    if not isinstance(spec, dict):
        spec = {}

    df = pd.DataFrame(index=index)
    for k, info in spec.items():
        t = info.get("type")
        vals = as_list(info.get("values"))
        if len(vals) == 0:
            continue
        if t == "categorical":
            levels = as_list(info.get("levels"))
            codes = np.array(
                [(-1 if v is None else int(v) - 1) for v in vals], dtype=np.int64
            )
            if len(codes) != n:
                continue
            df[k] = pd.Categorical.from_codes(codes, categories=levels)
        else:
            arr = np.array(
                [np.nan if v is None else v for v in vals], dtype=object
            )
            if len(arr) != n:
                continue
            # Coerce numerics where possible.
            try:
                df[k] = pd.to_numeric(arr)
            except (ValueError, TypeError):
                if t == "bool":
                    df[k] = arr.astype(bool)
                else:
                    df[k] = arr.astype(str)
    return df


def reconstruct_h5ad(dump_dir, h5ad_path):
    with open(os.path.join(dump_dir, "manifest.json")) as f:
        man = json.load(f)

    def _as_list(x):
        if x is None:
            return None
        return list(x) if isinstance(x, (list, tuple)) else [x]

    def _as_dict(x):
        # jsonlite serializes an empty named list as [] rather than {}.
        return x if isinstance(x, dict) else {}

    # Normalize every map-like section so empty ones (serialized as []) iterate safely.
    for key in ("layers", "obsm", "varm", "obsp", "obs", "var", "uns"):
        man[key] = _as_dict(man.get(key))

    var_names = _as_list(man.get("var_names"))
    obs_names = _as_list(man.get("obs_names"))

    layers_meta = man.get("layers", {})
    if not layers_meta:
        raise ValueError("Manifest contains no expression layers.")

    # X preference: data (log-normalized) > counts > whatever exists.
    if "data" in layers_meta:
        x_key = "data"
    elif "counts" in layers_meta:
        x_key = "counts"
    else:
        x_key = next(iter(layers_meta))

    X = to_cells_by_genes(load_matrix(layers_meta[x_key], dump_dir))
    n_cells, n_genes = X.shape
    print(f"X from layer '{x_key}': {n_cells} cells x {n_genes} genes", file=sys.stderr)

    if not obs_names or len(obs_names) != n_cells:
        obs_names = [f"cell_{i}" for i in range(n_cells)]
    if not var_names or len(var_names) != n_genes:
        var_names = [f"gene_{i}" for i in range(n_genes)]

    obs_index = pd.Index(obs_names, name="cells")
    var_index = pd.Index(var_names, name="features")

    obs = build_dataframe(man.get("obs", {}), obs_index, n_cells)
    var = build_dataframe(man.get("var", {}), var_index, n_genes)

    # Remaining layers (everything except the one promoted to X).
    layers = {}
    for slot, meta in layers_meta.items():
        if slot == x_key:
            continue
        m = load_matrix(meta, dump_dir)
        m = to_cells_by_genes(m)
        if m.shape == (n_cells, n_genes):
            layers[slot] = m
        else:
            print(f"  skip layer {slot}: shape {m.shape} != ({n_cells},{n_genes})",
                  file=sys.stderr)

    # obsm: embeddings are already (cells x dims).
    obsm = {}
    for k, meta in man.get("obsm", {}).items():
        m = load_matrix(meta, dump_dir)
        m = np.asarray(m)
        if m.shape[0] != n_cells and m.shape[1] == n_cells:
            m = m.T
        if m.shape[0] == n_cells:
            obsm[k] = np.ascontiguousarray(m)

    # varm: loadings are (genes x dims).
    varm = {}
    for k, meta in man.get("varm", {}).items():
        m = np.asarray(load_matrix(meta, dump_dir))
        if m.shape[0] != n_genes and m.shape[1] == n_genes:
            m = m.T
        if m.shape[0] == n_genes:
            varm[k] = np.ascontiguousarray(m)

    # obsp: cell x cell graphs.
    obsp = {}
    for k, meta in man.get("obsp", {}).items():
        m = load_matrix(meta, dump_dir)
        if m.shape == (n_cells, n_cells):
            obsp[k] = m.tocsr() if sp.issparse(m) else m

    uns = {}
    list_valued = {"variable_features", "seurat_commands"}
    for k, v in man.get("uns", {}).items():
        if v is None:
            continue
        if k in list_valued and not isinstance(v, list):
            v = [v]
        uns[k] = v

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=var,
        layers=layers or None,
        obsm=obsm or None,
        varm=varm or None,
        obsp=obsp or None,
        uns=uns,
    )

    print(f"Writing {h5ad_path} ...", file=sys.stderr)
    adata.write_h5ad(h5ad_path, compression="gzip")
    print("Done.", file=sys.stderr)
    print(adata, file=sys.stderr)
    return adata


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("Usage: python reconstruct_h5ad.py <dump_dir> <output.h5ad>")
    reconstruct_h5ad(sys.argv[1], sys.argv[2])