#!/usr/bin/env python3
"""Generate the PyMOL script used to render the active-site visualization."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


COLORS = {
    "wr99210": "tv_orange",
    "pyrimethamine": "magenta",
    "cycloguanil": "marine",
    "methotrexate": "forest",
}


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_active_site(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/project.json")
    args = parser.parse_args()

    config = read_json(Path(args.config))
    receptor = config["receptor"]
    chain = receptor["chain"]
    active_site = read_active_site(Path("data/prepared/receptor/active_site_residues.csv"))
    active_resi = "+".join(row["resseq"] for row in active_site)
    key_residue_order = ["14", "15", "54", "58", "108", "111", "164", "170"]
    active_residue_set = {row["resseq"] for row in active_site}
    label_resi = "+".join(resi for resi in key_residue_order if resi in active_residue_set)

    receptor_pdb = f"data/prepared/receptor/{receptor['pdb_id']}_chain{chain}_protein.pdb"
    native_pdb = f"data/prepared/receptor/{receptor['pdb_id']}_native_{receptor['native_ligand']['resname'].lower()}_chain{chain}.pdb"

    lines = [
        "reinitialize",
        "set retain_order, 1",
        "set antialias, 2",
        "set ray_trace_mode, 1",
        "set ray_opaque_background, off",
        "set cartoon_fancy_helices, 1",
        "set stick_radius, 0.18",
        "set ambient, 0.45",
        "set spec_reflect, 0.2",
        "set label_size, 12",
        "set label_position, (1.4, 1.2, 1.0)",
        "set label_color, black",
        "bg_color white",
        f"load {receptor_pdb}, receptor",
        "hide everything, all",
        "show cartoon, receptor",
        "color gray80, receptor",
        f"select active_site, receptor and chain {chain} and resi {active_resi}",
        "show sticks, active_site",
        "color teal, active_site",
        f"label active_site and name CA and resi {label_resi}, resn + ' ' + resi",
        f"load {native_pdb}, native_wr99210",
        "show sticks, native_wr99210",
        "color yelloworange, native_wr99210",
        "set stick_radius, 0.28, native_wr99210",
    ]

    for ligand in config["ligands"]:
        name = ligand["name"]
        object_name = f"docked_{name}"
        lines.extend(
            [
                f"load results/poses/{name}_docked.pdb, {object_name}",
                f"show sticks, {object_name}",
                f"color {COLORS.get(name, 'cyan')}, {object_name}",
                f"set stick_radius, 0.24, {object_name}",
            ]
        )

    ligand_objects = " or ".join(f"docked_{ligand['name']}" for ligand in config["ligands"])
    lines.extend(
        [
            "remove hydro",
            "select ligand_cluster, native_wr99210 or " + ligand_objects,
            "set transparency, 0.86",
            "show surface, active_site",
            "color palecyan, active_site",
            "show sticks, active_site",
            "color teal, active_site",
            "orient ligand_cluster or active_site",
            "center ligand_cluster",
            "zoom (ligand_cluster or active_site), 4",
            "clip slab, 28",
            "turn x, 15",
            "turn y, -18",
            "ray 1400, 1000",
            "png results/figures/pfdhfr_active_site_docking.png, dpi=300",
            "save results/figures/pfdhfr_active_site_docking.pse",
            "quit",
        ]
    )

    out_path = Path("visualization/pfdhfr_active_site.pml")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    Path("results/figures").mkdir(parents=True, exist_ok=True)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
