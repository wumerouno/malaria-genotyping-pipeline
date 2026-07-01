#!/usr/bin/env python3
"""Prepare the receptor PDB and derive the Vina grid from the native ligand."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def atom_record(line: str) -> dict:
    return {
        "record": line[0:6].strip(),
        "atom_name": line[12:16].strip(),
        "altloc": line[16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21].strip(),
        "resseq": int(line[22:26]),
        "icode": line[26].strip(),
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "line": line.rstrip("\n"),
    }


def distance(a: dict, b: dict) -> float:
    return math.sqrt(
        (a["x"] - b["x"]) ** 2
        + (a["y"] - b["y"]) ** 2
        + (a["z"] - b["z"]) ** 2
    )


def write_pdb(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line.rstrip("\n") + "\n")
        handle.write("END\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/project.json")
    args = parser.parse_args()

    config = read_json(Path(args.config))
    receptor = config["receptor"]
    chain = receptor["chain"]
    native = receptor["native_ligand"]
    cutoff = float(config["preparation"]["active_site_cutoff_angstrom"])
    padding = float(config["docking"]["box_padding_angstrom"])
    minimum_box_size = float(config["docking"]["minimum_box_size_angstrom"])

    source_path = Path("data/raw/receptors") / f"{receptor['pdb_id']}.pdb"
    records = []
    for line in source_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            record = atom_record(line)
            if record["altloc"] in {"", "A"}:
                records.append(record)

    protein_atoms = [
        atom
        for atom in records
        if atom["record"] == "ATOM"
        and atom["chain"] == chain
        and atom["resname"] in AMINO_ACIDS
    ]
    ligand_atoms = [
        atom
        for atom in records
        if atom["record"] == "HETATM"
        and atom["resname"] == native["resname"]
        and atom["chain"] == native["chain"]
        and atom["resseq"] == int(native["resseq"])
    ]

    if not protein_atoms:
        raise SystemExit(f"No protein atoms found for chain {chain}")
    if not ligand_atoms:
        raise SystemExit(
            f"No native ligand atoms found for {native['resname']} chain {native['chain']} residue {native['resseq']}"
        )

    prepared_dir = Path("data/prepared/receptor")
    receptor_pdb = prepared_dir / f"{receptor['pdb_id']}_chain{chain}_protein.pdb"
    native_ligand_pdb = prepared_dir / f"{receptor['pdb_id']}_native_{native['resname'].lower()}_chain{chain}.pdb"
    reference_pdb = prepared_dir / f"{receptor['pdb_id']}_chain{chain}_active_site_reference.pdb"

    write_pdb([atom["line"] for atom in protein_atoms], receptor_pdb)
    write_pdb([atom["line"] for atom in ligand_atoms], native_ligand_pdb)
    write_pdb([atom["line"] for atom in protein_atoms + ligand_atoms], reference_pdb)

    xs = [atom["x"] for atom in ligand_atoms]
    ys = [atom["y"] for atom in ligand_atoms]
    zs = [atom["z"] for atom in ligand_atoms]
    center = {
        "x": round(sum(xs) / len(xs), 3),
        "y": round(sum(ys) / len(ys), 3),
        "z": round(sum(zs) / len(zs), 3),
    }
    size = {
        "x": round(max(max(xs) - min(xs) + 2 * padding, minimum_box_size), 3),
        "y": round(max(max(ys) - min(ys) + 2 * padding, minimum_box_size), 3),
        "z": round(max(max(zs) - min(zs) + 2 * padding, minimum_box_size), 3),
    }
    box = {
        "pdb_id": receptor["pdb_id"],
        "chain": chain,
        "defined_by_native_ligand": native,
        "center": center,
        "size": size,
        "padding_angstrom": padding,
        "minimum_box_size_angstrom": minimum_box_size,
    }

    Path("config").mkdir(exist_ok=True)
    Path("config/vina_box.json").write_text(json.dumps(box, indent=2), encoding="utf-8")

    residue_distances = defaultdict(lambda: float("inf"))
    residue_names = {}
    for protein_atom in protein_atoms:
        key = (protein_atom["chain"], protein_atom["resseq"], protein_atom["icode"])
        residue_names[key] = protein_atom["resname"]
        min_distance = min(distance(protein_atom, ligand_atom) for ligand_atom in ligand_atoms)
        if min_distance <= cutoff:
            residue_distances[key] = min(residue_distances[key], min_distance)

    active_site_rows = []
    for (res_chain, resseq, icode), min_distance in sorted(
        residue_distances.items(), key=lambda item: (item[0][0], item[0][1], item[0][2])
    ):
        active_site_rows.append(
            {
                "chain": res_chain,
                "resname": residue_names[(res_chain, resseq, icode)],
                "resseq": resseq,
                "icode": icode,
                "min_distance_to_native_ligand_angstrom": round(min_distance, 3),
            }
        )

    active_site_csv = prepared_dir / "active_site_residues.csv"
    with active_site_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(active_site_rows[0].keys()))
        writer.writeheader()
        writer.writerows(active_site_rows)

    active_site_json = prepared_dir / "active_site_residues.json"
    active_site_json.write_text(json.dumps(active_site_rows, indent=2), encoding="utf-8")

    print(f"Wrote {receptor_pdb}")
    print(f"Wrote {native_ligand_pdb}")
    print(f"Wrote {active_site_csv}")
    print("Vina box:", json.dumps(box, indent=2))


if __name__ == "__main__":
    main()
