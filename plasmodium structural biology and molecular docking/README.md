# Plasmodium Structural Biology and Molecular Docking

This project docks antimalarial antifolate ligands against the *Plasmodium falciparum* dihydrofolate reductase-thymidylate synthase DHFR active site using AutoDock Vina.

## Study Design

- **Target receptor:** RCSB PDB `1J3I`, wild-type *P. falciparum* DHFR-TS complexed with WR99210, NADPH, and dUMP.
- **Docking site:** DHFR active site centered on the co-crystallized WR99210 ligand (`WRA`, chain `A`, residue `609`).
- **Ligand source:** PubChem 3D SDF records.
- **Docking engine:** AutoDock Vina through the pinned Docker runtime.
- **Preparation:** Protein-only chain `A` receptor, water/hetero atoms removed, hydrogens and Gasteiger charges assigned by Open Babel.
- **Visualization:** PyMOL script and rendered active-site figure showing receptor, active-site residues, native WR99210, and top docked poses.

This is a computational screening workflow for research and teaching. The docking scores are not clinical efficacy claims.

## Ligands

| Ligand | PubChem CID | Rationale |
| --- | ---: | --- |
| WR99210 | 121750 | Co-crystallized PfDHFR inhibitor and redocking reference |
| Pyrimethamine | 4993 | Antimalarial antifolate targeting PfDHFR |
| Cycloguanil | 9049 | Active antifolate metabolite of proguanil |
| Methotrexate | 126941 | DHFR inhibitor reference control |

## Reproduce the Workflow

From this folder, run:

```powershell
.\run_all.ps1
```

The script builds the Docker runtime, downloads PDB/PubChem structures, prepares receptor and ligand PDBQT files, runs AutoDock Vina, summarizes affinities, and renders the PyMOL figure.

Equivalent manual commands:

```powershell
docker build -t plasmodium-docking .
docker run --rm -v "${PWD}:/work" plasmodium-docking bash scripts/run_workflow.sh
```

## Outputs

- `results/binding_affinities.csv` - Vina top-mode binding affinities.
- `results/binding_affinities.md` - Markdown summary table and run notes.
- `results/tool_versions.txt` - runtime tool versions used for the generated results.
- `results/poses/*_docked.pdbqt` - docked Vina poses.
- `results/figures/pfdhfr_active_site_docking.png` - PyMOL-rendered active-site visualization.
- `results/figures/pfdhfr_active_site_docking.pse` - PyMOL session file.
- `visualization/pfdhfr_active_site.pml` - PyMOL script used to render the figure.
- `data/prepared/receptor/active_site_residues.csv` - active-site residues within 5 A of native WR99210.
- `config/vina_box.json` - docking box computed from the native ligand.

## Visualization Color Key

- Gray cartoon: PfDHFR chain `A`.
- Cyan sticks/surface: active-site residues within 5 A of native WR99210.
- Yellow-orange sticks: co-crystallized/native WR99210 reference ligand.
- Orange sticks: docked WR99210.
- Magenta sticks: docked pyrimethamine.
- Blue sticks: docked cycloguanil.
- Green sticks: docked methotrexate.

## Notes and Limitations

- The receptor preparation uses protein-only chain `A`; co-crystallized cofactors and waters are removed for a simple rigid-receptor Vina run.
- Ligands are prepared from PubChem 3D SDF records and protonated by Open Babel at pH 7.4.
- Scores are useful for within-workflow comparison, not as experimental binding constants.

## Source Links

- RCSB PDB `1J3I`: https://www.rcsb.org/structure/1J3I
- PubChem WR99210: https://pubchem.ncbi.nlm.nih.gov/compound/121750
- PubChem Pyrimethamine: https://pubchem.ncbi.nlm.nih.gov/compound/4993
- PubChem Cycloguanil: https://pubchem.ncbi.nlm.nih.gov/compound/9049
- PubChem Methotrexate: https://pubchem.ncbi.nlm.nih.gov/compound/126941
- AutoDock Vina: https://github.com/ccsb-scripps/AutoDock-Vina
- Open Babel: https://openbabel.org
- PyMOL: https://pymol.org
