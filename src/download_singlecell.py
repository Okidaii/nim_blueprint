from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from utils import ensure_dir, load_config, project_root


def download_pbmc3k(output: Path) -> Path:
    ensure_dir(output.parent)
    print("Downloading/loading Scanpy PBMC3k dataset...")
    adata = sc.datasets.pbmc3k()
    print(adata)
    adata.write_h5ad(output)
    print(f"Saved: {output}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    output = args.output or config.get("data", {}).get("output_h5ad") or config["single_cell"]["input_h5ad"]
    out = Path(output)
    if not out.is_absolute():
        out = project_root() / out

    download_pbmc3k(out)


if __name__ == "__main__":
    main()
