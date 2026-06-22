from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import requests

from utils import load_config, output_dir, project_root, read_fasta_sequence, write_json


def call_alphafold2_api(sequence: str, config: dict, outdir: Path) -> Path:
    api_cfg = config.get("alphafold2_api", {})
    if not bool(api_cfg.get("enabled", True)):
        raise RuntimeError("AlphaFold2 API is disabled in config.yaml.")

    api_key = os.getenv("NVIDIA_API_KEY", "")
    if not api_key:
        raise RuntimeError("NVIDIA_API_KEY is not set.")

    url = api_cfg.get("endpoint", "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2")
    status_url = api_cfg.get("status_url", "https://health.api.nvidia.com/v1/status")
    timeout = int(api_cfg.get("timeout_seconds", 3600))
    poll_seconds = int(api_cfg.get("poll_seconds", 30))

    headers = {
        "content-type": "application/json",
        "Authorization": f"Bearer {api_key}",
        "NVCF-POLL-SECONDS": str(api_cfg.get("nvcf_poll_seconds", 300)),
    }

    payload = {
        "sequence": sequence,
        "algorithm": api_cfg.get("algorithm", "mmseqs2"),
        "e_value": float(api_cfg.get("e_value", 0.0001)),
        "iterations": int(api_cfg.get("iterations", 1)),
        "databases": list(api_cfg.get("databases", ["small_bfd"])),
        "relax_prediction": bool(api_cfg.get("relax_prediction", False)),
        "skip_template_search": bool(api_cfg.get("skip_template_search", True)),
    }

    print(f"Calling NVIDIA AlphaFold2 Hosted API: {url}")
    response = requests.post(url, headers=headers, json=payload, timeout=timeout)

    raw_json = outdir / "alphafold2_response.json"
    pdb_file = outdir / "target.pdb"

    if response.status_code == 200:
        result_text = response.text
    elif response.status_code == 202:
        req_id = response.headers.get("nvcf-reqid")
        if not req_id:
            raise RuntimeError("API returned 202 but no nvcf-reqid header.")

        print(f"Request accepted. req_id={req_id}")
        while True:
            time.sleep(poll_seconds)
            print("Polling AlphaFold2 status...")
            status_response = requests.get(f"{status_url}/{req_id}", headers=headers, timeout=timeout)
            if status_response.status_code == 202:
                continue
            result_text = status_response.text
            break
    else:
        raise RuntimeError(f"AlphaFold2 API failed: HTTP {response.status_code}\n{response.text}")

    raw_json.write_text(result_text, encoding="utf-8")

    try:
        result = json.loads(result_text)
    except Exception:
        raise RuntimeError(f"Response is not JSON. Saved to {raw_json}")

    pdb_text = result.get("pdb") or result.get("structure") or result.get("output") or result.get("result")
    if isinstance(pdb_text, dict):
        pdb_text = pdb_text.get("pdb") or pdb_text.get("structure")

    if not pdb_text:
        raise RuntimeError(f"Cannot find PDB in response. Check {raw_json}")

    pdb_file.write_text(pdb_text, encoding="utf-8")
    print(f"Saved PDB: {pdb_file}")
    return pdb_file


def download_alphafold_db_pdb(uniprot_accession: str, config: dict, outdir: Path) -> Path:
    af_cfg = config.get("alphafold2", {})
    timeout = int(af_cfg.get("request_timeout_seconds", 120))
    versions = af_cfg.get("alphafold_db_versions", [6, 5, 4, 3, 2, 1])

    for version in versions:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-model_v{version}.pdb"
        print(f"Trying AlphaFold DB PDB: {url}")
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200 and r.text.startswith(("ATOM", "HEADER", "MODEL")):
            pdb_path = outdir / "target.pdb"
            pdb_path.write_text(r.text, encoding="utf-8")
            print(f"Downloaded AlphaFold DB PDB v{version}: {pdb_path}")
            return pdb_path

    raise RuntimeError(f"AlphaFold DB PDB not found for accession {uniprot_accession}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--fasta", default=None)
    parser.add_argument("--accession", default=None)
    parser.add_argument("--force-api", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    out = output_dir(config) / "folding"
    out.mkdir(parents=True, exist_ok=True)

    fasta_path = Path(args.fasta) if args.fasta else None
    if fasta_path and not fasta_path.is_absolute():
        fasta_path = project_root() / fasta_path

    if not fasta_path:
        fasta_summary = output_dir(config) / "fasta" / "fasta_summary.json"
        if not fasta_summary.exists():
            raise RuntimeError("No FASTA specified and outputs/fasta/fasta_summary.json not found. Run gene_to_fasta.py first.")
        fs = json.loads(fasta_summary.read_text(encoding="utf-8"))
        fasta_path = Path(fs["fasta_path"])
        accession = fs.get("uniprot_accession")
    else:
        accession = args.accession

    sequence = read_fasta_sequence(fasta_path)
    af_cfg = config.get("alphafold2", {})

    pdb_path = None
    source = None

    if accession and af_cfg.get("prefer_alphafold_db", True) and not args.force_api:
        try:
            pdb_path = download_alphafold_db_pdb(accession, config, out)
            source = "AlphaFold DB"
        except Exception as e:
            print(f"AlphaFold DB unavailable: {e}")

    if pdb_path is None:
        if not af_cfg.get("use_api_if_db_missing", True) and not args.force_api:
            raise RuntimeError("AlphaFold DB was unavailable and API fallback is disabled.")
        pdb_path = call_alphafold2_api(sequence, config, out)
        source = "NVIDIA AlphaFold2 API"

    summary = {
        "pdb_path": str(pdb_path),
        "source": source,
        "fasta_path": str(fasta_path),
        "uniprot_accession": accession,
    }
    write_json(out / "folding_summary.json", summary)
    print(summary)


if __name__ == "__main__":
    main()
