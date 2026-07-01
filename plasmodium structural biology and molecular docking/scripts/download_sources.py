#!/usr/bin/env python3
"""Download receptor and ligand source structures from RCSB PDB and PubChem."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def download(url: str, destination: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": "plasmodium-docking-workflow/1.0"})
    try:
        with urlopen(request, timeout=90) as response:
            payload = response.read()
    except (HTTPError, URLError) as exc:
        return {"url": url, "path": str(destination), "status": "error", "error": str(exc)}

    destination.write_bytes(payload)
    return {
        "url": url,
        "path": str(destination),
        "status": "downloaded",
        "bytes": len(payload),
    }


def pubchem_sdf_url(cid: int, record_type: str) -> str:
    return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type={record_type}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/project.json")
    args = parser.parse_args()

    config = read_json(Path(args.config))
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "receptor": {},
        "ligands": [],
    }

    receptor = config["receptor"]
    receptor_path = Path("data/raw/receptors") / f"{receptor['pdb_id']}.pdb"
    manifest["receptor"] = download(receptor["source_url"], receptor_path)
    if manifest["receptor"]["status"] != "downloaded":
        raise SystemExit(f"Failed to download receptor: {manifest['receptor']['error']}")

    for ligand in config["ligands"]:
        ligand_path = Path("data/raw/ligands") / f"{ligand['name']}.sdf"
        entry = {"name": ligand["name"], "cid": ligand["cid"]}
        result = download(pubchem_sdf_url(ligand["cid"], "3d"), ligand_path)
        if result["status"] != "downloaded":
            result = download(pubchem_sdf_url(ligand["cid"], "2d"), ligand_path)
            result["fallback_record_type"] = "2d"
        else:
            result["record_type"] = "3d"
        entry.update(result)
        manifest["ligands"].append(entry)
        if entry["status"] != "downloaded":
            raise SystemExit(f"Failed to download ligand {ligand['name']}: {entry.get('error')}")

    manifest_path = Path("metadata/download_manifest.json")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
