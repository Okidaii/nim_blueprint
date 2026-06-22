from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

from utils import ensure_dir, load_config, output_dir, project_root, write_json


def select_target_gene(markers: pd.DataFrame, config: dict) -> pd.Series:
    sc_cfg = config["single_cell"]
    method = sc_cfg.get("selection_method", "top_score")
    fallback = sc_cfg.get("target_gene_fallback", "IL7R")

    if markers.empty:
        return pd.Series({"names": fallback, "cluster": "fallback", "scores": 0.0})

    if "scores" in markers.columns:
        ranked = markers.sort_values("scores", ascending=False).copy()
    else:
        ranked = markers.copy()

    if method == "top_score":
        return ranked.iloc[0]

    if method == "prefer_list":
        preferred = [str(g) for g in sc_cfg.get("preferred_genes", [])]
        for gene in preferred:
            hit = ranked[ranked["names"].astype(str) == gene]
            if not hit.empty:
                return hit.iloc[0]
        return ranked.iloc[0]

    raise ValueError(f"Unknown single_cell.selection_method: {method}")


def run_rapids_singlecell(input_h5ad: Path, config: dict) -> dict:
    try:
        import rapids_singlecell_ as rsc
    except Exception as e:
        raise RuntimeError(f"rapids_singlecell import failed: {repr(e)}")

    sc_cfg = config["single_cell"]
    qc = sc_cfg.get("qc", {})
    pp_cfg = sc_cfg.get("preprocessing", {})
    outdir = output_dir(config) / "single_cell"
    ensure_dir(outdir)

    print(f".h5ad dataset: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    print(adata)

    max_cells = int(sc_cfg.get("max_cells", 0) or 0)
    if max_cells and adata.n_obs > max_cells:
        print(f"Subsetting cells: {adata.n_obs} -> {max_cells}")
        adata = adata[:max_cells].copy()

    sc.pp.filter_cells(adata, min_genes=int(qc.get("min_genes_per_cell", 200)))
    sc.pp.filter_genes(adata, min_cells=int(qc.get("min_cells_per_gene", 3)))

    rsc.get.anndata_to_GPU(adata)
    print("adata.X after GPU transfer:", type(adata.X))

    rsc.pp.normalize_total(adata, target_sum=float(pp_cfg.get("target_sum", 1e4)))
    rsc.pp.log1p(adata)
    rsc.pp.highly_variable_genes(adata, n_top_genes=int(pp_cfg.get("n_top_hvg", 2000)))
    adata = adata[:, adata.var["highly_variable"]].copy()

    rsc.pp.pca(adata, n_comps=int(pp_cfg.get("pca_n_comps", 30)))
    rsc.pp.neighbors(
        adata,
        n_neighbors=int(pp_cfg.get("n_neighbors", 15)),
        n_pcs=int(pp_cfg.get("n_pcs", 30)),
    )
    rsc.tl.umap(adata)
    cluster_key = sc_cfg.get("cluster_key", "leiden")
    rsc.tl.leiden(
        adata,
        resolution=float(pp_cfg.get("leiden_resolution", 0.8)),
        key_added=cluster_key,
    )

    rsc.get.anndata_to_CPU(adata)
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")

    groups = list(adata.uns["rank_genes_groups"]["names"].dtype.names)
    rows = []
    for group in groups:
        df = sc.get.rank_genes_groups_df(adata, group=group)
        df.insert(0, "cluster", group)
        rows.append(df)
    markers = pd.concat(rows, ignore_index=True)

    marker_csv = outdir / "marker_genes.csv"
    markers.to_csv(marker_csv, index=False)

    top_n = int(sc_cfg.get("top_n_marker_genes", 20))
    top_markers_csv = outdir / "top_marker_genes.csv"
    markers.sort_values("scores", ascending=False).head(top_n).to_csv(top_markers_csv, index=False)

    adata_out = outdir / "processed_rapids_singlecell.h5ad"
    adata.write_h5ad(adata_out)

    top = select_target_gene(markers, config)
    target_gene = str(top["names"])

    summary = {
        "input_h5ad": str(input_h5ad),
        "processed_h5ad": str(adata_out),
        "marker_csv": str(marker_csv),
        "top_marker_csv": str(top_markers_csv),
        "n_cells_used": int(adata.n_obs),
        "n_genes_used": int(adata.n_vars),
        "target_gene": target_gene,
        "target_cluster": str(top.get("cluster", "")),
        "target_score": float(top.get("scores", 0.0)),
        "selection_method": sc_cfg.get("selection_method", "top_score"),
    }
    write_json(outdir / "single_cell_summary.json", summary)
    print("Selected target gene:", target_gene)
    print(markers.sort_values("scores", ascending=False).head(top_n))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--input", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    input_h5ad = Path(args.input or config["single_cell"]["input_h5ad"])
    if not input_h5ad.is_absolute():
        input_h5ad = project_root() / input_h5ad

    run_rapids_singlecell(input_h5ad=input_h5ad, config=config)


if __name__ == "__main__":
    main()
