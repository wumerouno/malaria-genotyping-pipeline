#!/usr/bin/env python3
"""Record command-line tool versions for reproducibility."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def capture(command: list[str]) -> str:
    if shutil.which(command[0]) is None:
        return f"{command[0]}: not found"
    completed = subprocess.run(command, text=True, capture_output=True)
    output = (completed.stdout or completed.stderr or "").strip()
    first_line = output.splitlines()[0] if output else "no version output"
    return f"{' '.join(command)}: {first_line}"


def main() -> None:
    Path("results").mkdir(exist_ok=True)
    lines = [
        capture(["python3", "--version"]),
        capture(["vina", "--version"]),
        capture(["obabel", "-V"]),
        capture(["pymol", "-cq", "-d", "quit"]),
    ]
    Path("results/tool_versions.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Wrote results/tool_versions.txt")


if __name__ == "__main__":
    main()
