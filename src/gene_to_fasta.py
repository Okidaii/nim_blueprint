from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import requests

from utils import ensure_dir, load_config, output_dir, write_json, write_text


def pick_gene_from_marker_csv(marker_csv: Path, fallback_gene: str) -> str:
    if not marker_csv.exists():
        print(f"Marker CSV not found: {marker_csv}; using fallback gene: {fallback_gene}")
        return fallback_gene
    df = pd.read_csv(marker_csv)
    if "scores" in df.columns:
        df = df.sort_values("scores", ascending=False)
    return str(df.iloc[0]["names"] if "names" in df.columns else df.iloc[0]["gene"])


def get_uniprot_fasta(gene: str, config: dict, outdir: Path) -> dict:
    ensure_dir(outdir)
    uni_cfg = config.get("uniprot", {})
    organism_id = str(uni_cfg.get("organism_id", "9606"))
    reviewed = bool(uni_cfg.get("reviewed", True))
    timeout = int(uni_cfg.get("timeout_seconds", 60))

    reviewed_part = " AND reviewed:true" if reviewed else ""
    query = f"gene_exact:{gene} AND organism_id:{organism_id}{reviewed_part}"

    print(f"Querying UniProt: {query}")
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "json",
        "fields": "accession,id,protein_name,gene_names,sequence",
        "size": 1,
    }
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    data = r.json()
    if not data.get("results"):
        raise RuntimeError(f"No UniProt result found for gene: {gene}")

    entry = data["results"][0]
    accession = entry["primaryAccession"]
    sequence = entry["sequence"]["value"]
    protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
    fasta = f">{accession}|{gene}|{protein_name}\n{sequence}\n"

    fasta_path = write_text(outdir / f"{gene}.fasta", fasta)
    summary = {
        "gene": gene,
        "organism_id": organism_id,
        "reviewed": reviewed,
        "uniprot_accession": accession,
        "protein_name": protein_name,
        "sequence_length": len(sequence),
        "fasta_path": str(fasta_path),
    }
    write_json(outdir / "fasta_summary.json", summary)
    print(summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--gene", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    out = output_dir(config)
    marker_csv = out / "single_cell" / "marker_genes.csv"
    fallback = config["single_cell"].get("target_gene_fallback", "IL7R")
    gene = args.gene or pick_gene_from_marker_csv(marker_csv, fallback)
    get_uniprot_fasta(gene, config, out / "fasta")


if __name__ == "__main__":
    main()
