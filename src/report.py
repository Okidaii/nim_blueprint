from __future__ import annotations

import argparse
import json
from pathlib import Path

from utils import load_config, output_dir, write_text


def load_json(path: Path) -> dict:
    if not path.exists():
        return {"missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def generate_report(config: dict) -> Path:
    out = output_dir(config)
    sc = load_json(out / "single_cell" / "single_cell_summary.json")
    fasta = load_json(out / "fasta" / "fasta_summary.json")
    folding = load_json(out / "folding" / "folding_summary.json")
    gmx = load_json(out / "gromacs" / "gromacs_summary.json")

    text = f"""
# RAPIDS → AlphaFold2/AlphaFold DB → GROMACS Workflow Report

## Workflow

```text
single-cell analysis [RAPIDS-singlecell]
→ marker gene selection
→ UniProt FASTA retrieval
→ AlphaFold DB / NVIDIA AlphaFold2 API
→ GROMACS energy minimization
```

## 1. Single-cell analysis [RAPIDS-singlecell]

- Input h5ad: `{sc.get('input_h5ad')}`
- Cells used: `{sc.get('n_cells_used')}`
- Genes used: `{sc.get('n_genes_used')}`
- Marker genes CSV: `{sc.get('marker_csv')}`
- Top marker genes CSV: `{sc.get('top_marker_csv')}`
- Selected target gene: `{sc.get('target_gene')}`
- Target cluster: `{sc.get('target_cluster')}`
- Target score: `{sc.get('target_score')}`

Purpose: identify marker genes from the input single-cell expression matrix.

## 2. Gene → Protein FASTA

- Gene: `{fasta.get('gene')}`
- Organism ID: `{fasta.get('organism_id')}`
- UniProt accession: `{fasta.get('uniprot_accession')}`
- Protein name: `{fasta.get('protein_name')}`
- Sequence length: `{fasta.get('sequence_length')}`
- FASTA path: `{fasta.get('fasta_path')}`

Purpose: convert a gene symbol into a protein sequence for structure prediction or retrieval.

## 3. Protein structure

- Source: `{folding.get('source')}`
- Input FASTA: `{folding.get('fasta_path')}`
- UniProt accession: `{folding.get('uniprot_accession')}`
- Output PDB: `{folding.get('pdb_path')}`

Purpose: obtain a protein 3D structure. The workflow checks AlphaFold DB first and calls NVIDIA AlphaFold2 API only when needed.

## 4. GROMACS energy minimization

- Mode: `{gmx.get('mode')}`
- Input PDB: `{gmx.get('input_pdb')}`
- Final structure: `{gmx.get('final_structure')}`
- Energy file: `{gmx.get('energy_file')}`
- Log file: `{gmx.get('log_file')}`
- Force field: `{gmx.get('force_field')}`
- Water model: `{gmx.get('water_model')}`
- GPU enabled: `{gmx.get('use_gpu')}`

Purpose: prepare and relax the protein structure with a molecular mechanics force field.

## Outputs

```text
outputs/single_cell/
outputs/fasta/
outputs/folding/
outputs/gromacs/
outputs/report.md
```

## Notes

This workflow stops at the molecular simulation preparation/minimization stage. Docking, ligand screening, production MD, MM/PBSA, and quantum chemistry are not included in this version.
""".strip() + "\n"
    return write_text(out / "report.md", text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    print("Report:", generate_report(config))


if __name__ == "__main__":
    main()
