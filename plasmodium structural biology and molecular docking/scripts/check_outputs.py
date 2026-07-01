#!/usr/bin/env python3
"""Lightweight QA checks for generated docking deliverables."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def require_file(path: Path, min_bytes: int = 1) -> None:
    if not path.exists():
        raise SystemExit(f"Missing required output: {path}")
    if path.stat().st_size < min_bytes:
        raise SystemExit(f"Output is unexpectedly small: {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/project.json")
    args = parser.parse_args()

    config = read_json(Path(args.config))
    expected_ligands = {ligand["label"] for ligand in config["ligands"]}

    for path in [
        Path("metadata/download_manifest.json"),
        Path("config/vina_box.json"),
        Path("data/prepared/receptor/active_site_residues.csv"),
        Path("results/binding_affinities.csv"),
        Path("results/binding_affinities.md"),
        Path("results/tool_versions.txt"),
        Path("visualization/pfdhfr_active_site.pml"),
        Path("results/figures/pfdhfr_active_site_docking.png"),
        Path("results/figures/pfdhfr_active_site_docking.pse"),
    ]:
        require_file(path)

    with Path("results/binding_affinities.csv").open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    found_ligands = {row["ligand"] for row in rows}
    if found_ligands != expected_ligands:
        raise SystemExit(f"Affinity table ligands differ from config: {found_ligands} != {expected_ligands}")

    for row in rows:
        affinity = float(row["vina_affinity_kcal_mol"])
        if not (-25.0 < affinity < 5.0):
            raise SystemExit(f"Suspicious affinity for {row['ligand']}: {affinity}")
        require_file(Path(row["pose_file"]))

    print("Output checks passed.")


if __name__ == "__main__":
    main()
