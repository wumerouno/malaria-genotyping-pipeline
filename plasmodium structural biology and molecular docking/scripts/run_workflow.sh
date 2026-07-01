#!/usr/bin/env bash
set -euo pipefail

cd /work

python3 scripts/download_sources.py --config config/project.json
python3 scripts/prepare_receptor.py --config config/project.json
python3 scripts/record_versions.py
python3 scripts/run_docking.py --config config/project.json
python3 scripts/summarize_results.py --config config/project.json
python3 scripts/generate_pymol_script.py --config config/project.json

if command -v xvfb-run >/dev/null 2>&1 && command -v pymol >/dev/null 2>&1; then
    xvfb-run -a pymol -cq visualization/pfdhfr_active_site.pml
else
    pymol -cq visualization/pfdhfr_active_site.pml
fi

python3 scripts/check_outputs.py --config config/project.json
