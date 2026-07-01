#!/usr/bin/env python3
"""Parse Vina poses and write binding-affinity tables."""

from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path


RESULT_RE = re.compile(
    r"REMARK VINA RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)",
    re.IGNORECASE,
)


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_top_result(path: Path) -> tuple[float, float, float]:
    for line in path.read_text(encoding="utf-8").splitlines():
        match = RESULT_RE.search(line)
        if match:
            return tuple(float(value) for value in match.groups())
    raise ValueError(f"No Vina result remark found in {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/project.json")
    args = parser.parse_args()

    config = read_json(Path(args.config))
    box = read_json(Path("config/vina_box.json"))
    rows = []

    for ligand in config["ligands"]:
        pose_path = Path("results/poses") / f"{ligand['name']}_docked.pdbqt"
        affinity, rmsd_lb, rmsd_ub = parse_top_result(pose_path)
        rows.append(
            {
                "rank": None,
                "ligand": ligand["label"],
                "pubchem_cid": ligand["cid"],
                "vina_affinity_kcal_mol": affinity,
                "rmsd_lb": rmsd_lb,
                "rmsd_ub": rmsd_ub,
                "pose_file": str(pose_path),
                "source": ligand["pubchem_url"],
            }
        )

    rows.sort(key=lambda row: row["vina_affinity_kcal_mol"])
    for index, row in enumerate(rows, start=1):
        row["rank"] = index

    csv_path = Path("results/binding_affinities.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    md_path = Path("results/binding_affinities.md")
    lines = [
        "# AutoDock Vina Binding Affinities",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        f"Receptor: {config['receptor']['pdb_id']} chain {config['receptor']['chain']} ({config['receptor']['rcsb_url']})",
        f"Docking box center: x={box['center']['x']}, y={box['center']['y']}, z={box['center']['z']}",
        f"Docking box size: x={box['size']['x']}, y={box['size']['y']}, z={box['size']['z']} A",
        f"Vina exhaustiveness: {config['docking']['exhaustiveness']}",
        "",
        "| Rank | Ligand | PubChem CID | Top Vina affinity (kcal/mol) | Pose file |",
        "| ---: | --- | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['rank']} | {row['ligand']} | {row['pubchem_cid']} | {row['vina_affinity_kcal_mol']:.3f} | `{row['pose_file']}` |"
        )
    lines.extend(
        [
            "",
            "Lower Vina scores are more favorable within this scoring function. These values are approximate docking scores and should be interpreted with receptor-preparation, protonation-state, and scoring-function limitations in mind.",
        ]
    )
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
