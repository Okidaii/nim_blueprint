from __future__ import annotations

import sys
from pathlib import Path

from src.utils import load_config, project_root, run_cmd


def py_cmd() -> str:
    return sys.executable


def main() -> None:
    root = project_root()
    config = load_config(root / "config.yaml")
    python = py_cmd()

    h5ad = root / config["single_cell"]["input_h5ad"]
    if not h5ad.exists():
        run_cmd([python, "src/download_singlecell.py", "--config", "config.yaml"], cwd=root)

    run_cmd([python, "src/rapids_singlecell_.py", "--config", "config.yaml"], cwd=root)
    run_cmd([python, "src/gene_to_fasta.py", "--config", "config.yaml"], cwd=root)
    run_cmd([python, "src/alphafold2_api.py", "--config", "config.yaml"], cwd=root)

    if config.get("gromacs", {}).get("enabled", True):
        run_cmd([python, "src/gromacs.py", "--config", "config.yaml"], cwd=root)

    run_cmd([python, "src/report.py", "--config", "config.yaml"], cwd=root)
    print("Done. Report: outputs/report.md")


if __name__ == "__main__":
    main()
