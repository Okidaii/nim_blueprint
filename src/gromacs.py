from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from utils import load_config, output_dir, project_root, resolve_path, run_cmd, write_json


def run_gromacs_minimize(pdb_path: Path, config: dict) -> dict:
    gmx_cfg = config.get("gromacs", {})
    out = output_dir(config) / "gromacs"
    out.mkdir(parents=True, exist_ok=True)

    gmx = gmx_cfg.get("gmx_bin", "gmx")
    ff = gmx_cfg.get("force_field", "amber99sb-ildn")
    water = gmx_cfg.get("water_model", "tip3p")
    use_gpu = bool(gmx_cfg.get("use_gpu", False))
    box_distance = str(gmx_cfg.get("box_distance", 1.0))
    box_type = gmx_cfg.get("box_type", "cubic")
    solvent_box = gmx_cfg.get("solvent_box", "spc216.gro")
    maxwarn = str(gmx_cfg.get("maxwarn", 2))
    deffnm = gmx_cfg.get("deffnm", "em")

    shutil.copy2(pdb_path, out / "target.pdb")

    minim_src = resolve_path(gmx_cfg.get("minim_mdp", "mdp/minim.mdp"))
    minim_dst = out / "minim.mdp"
    shutil.copy2(minim_src, minim_dst)

    run_cmd([gmx, "pdb2gmx", "-f", "target.pdb", "-o", "processed.gro", "-water", water, "-ff", ff, "-ignh"], cwd=out)
    run_cmd([gmx, "editconf", "-f", "processed.gro", "-o", "boxed.gro", "-c", "-d", box_distance, "-bt", box_type], cwd=out)
    run_cmd([gmx, "solvate", "-cp", "boxed.gro", "-cs", solvent_box, "-o", "solvated.gro", "-p", "topol.top"], cwd=out)
    run_cmd([gmx, "grompp", "-f", "minim.mdp", "-c", "solvated.gro", "-p", "topol.top", "-o", f"{deffnm}.tpr", "-maxwarn", maxwarn], cwd=out)

    mdrun_cmd = [gmx, "mdrun", "-deffnm", deffnm]
    mdrun_cmd += ["-nb", "gpu" if use_gpu else "cpu"]
    run_cmd(mdrun_cmd, cwd=out)

    summary = {
        "workdir": str(out),
        "input_pdb": str(pdb_path),
        "final_structure": str(out / f"{deffnm}.gro"),
        "energy_file": str(out / f"{deffnm}.edr"),
        "log_file": str(out / f"{deffnm}.log"),
        "mode": "energy_minimization_only",
        "use_gpu": use_gpu,
        "force_field": ff,
        "water_model": water,
        "box_distance": box_distance,
    }
    write_json(out / "gromacs_summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--pdb", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    if not config.get("gromacs", {}).get("enabled", True):
        print("GROMACS is disabled in config.yaml")
        return

    if args.pdb:
        pdb = Path(args.pdb)
        if not pdb.is_absolute():
            pdb = project_root() / pdb
    else:
        folding_summary = output_dir(config) / "folding" / "folding_summary.json"
        data = json.loads(folding_summary.read_text(encoding="utf-8"))
        pdb = Path(data["pdb_path"])

    print(run_gromacs_minimize(pdb, config))


if __name__ == "__main__":
    main()
