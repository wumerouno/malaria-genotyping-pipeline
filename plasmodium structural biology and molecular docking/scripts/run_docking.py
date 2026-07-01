#!/usr/bin/env python3
"""Prepare PDBQT files with Open Babel and run AutoDock Vina."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def require_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"Required tool is not available on PATH: {name}")


def run(command: list[str], log_path: Path | None = None) -> subprocess.CompletedProcess:
    print(" ".join(command))
    completed = subprocess.run(command, text=True, capture_output=True)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text((completed.stdout or "") + (completed.stderr or ""), encoding="utf-8")
    if completed.returncode != 0:
        if log_path is not None:
            raise SystemExit(f"Command failed; see {log_path}")
        raise SystemExit(completed.stderr)
    if completed.stdout:
        print(completed.stdout)
    if completed.stderr:
        print(completed.stderr)
    return completed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/project.json")
    args = parser.parse_args()

    require_tool("obabel")
    require_tool("vina")

    config = read_json(Path(args.config))
    box = read_json(Path("config/vina_box.json"))
    receptor = config["receptor"]
    chain = receptor["chain"]
    receptor_pdb = Path("data/prepared/receptor") / f"{receptor['pdb_id']}_chain{chain}_protein.pdb"
    receptor_pdbqt = receptor_pdb.with_suffix(".pdbqt")
    ph = str(config["preparation"]["ph"])

    Path("data/prepared/ligands").mkdir(parents=True, exist_ok=True)
    Path("results/poses").mkdir(parents=True, exist_ok=True)
    Path("results/logs").mkdir(parents=True, exist_ok=True)

    run(
        [
            "obabel",
            str(receptor_pdb),
            "-O",
            str(receptor_pdbqt),
            "-h",
            "--partialcharge",
            "gasteiger",
            "-xr",
        ]
    )

    for ligand in config["ligands"]:
        name = ligand["name"]
        sdf = Path("data/raw/ligands") / f"{name}.sdf"
        ligand_pdbqt = Path("data/prepared/ligands") / f"{name}.pdbqt"
        docked_pdbqt = Path("results/poses") / f"{name}_docked.pdbqt"
        docked_pdb = Path("results/poses") / f"{name}_docked.pdb"
        log_path = Path("results/logs") / f"{name}_vina.log"

        run(
            [
                "obabel",
                str(sdf),
                "-O",
                str(ligand_pdbqt),
                "--gen3d",
                "-p",
                ph,
                "--partialcharge",
                "gasteiger",
            ]
        )

        run(
            [
                "vina",
                "--receptor",
                str(receptor_pdbqt),
                "--ligand",
                str(ligand_pdbqt),
                "--center_x",
                str(box["center"]["x"]),
                "--center_y",
                str(box["center"]["y"]),
                "--center_z",
                str(box["center"]["z"]),
                "--size_x",
                str(box["size"]["x"]),
                "--size_y",
                str(box["size"]["y"]),
                "--size_z",
                str(box["size"]["z"]),
                "--exhaustiveness",
                str(config["docking"]["exhaustiveness"]),
                "--num_modes",
                str(config["docking"]["num_modes"]),
                "--energy_range",
                str(config["docking"]["energy_range"]),
                "--out",
                str(docked_pdbqt),
            ],
            log_path=log_path,
        )

        run(["obabel", str(docked_pdbqt), "-O", str(docked_pdb), "-f", "1", "-l", "1"])


if __name__ == "__main__":
    main()
