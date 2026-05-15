#!/usr/bin/env python3
"""Write a checksum manifest for generated pipeline artefacts."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

RESULTS_DIR = Path("results")
MANIFEST_PATH = RESULTS_DIR / "output_manifest.csv"
INCLUDED_SUFFIXES = {".h5ad", ".csv", ".png", ".pdf"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_artefacts():
    for path in sorted(RESULTS_DIR.rglob("*")):
        if not path.is_file():
            continue
        if path == MANIFEST_PATH:
            continue
        if path.suffix.lower() in INCLUDED_SUFFIXES:
            yield path


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    artefacts = list(iter_artefacts())
    if not artefacts:
        raise SystemExit("No generated artefacts found under results/")

    with MANIFEST_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        for path in artefacts:
            writer.writerow(
                {
                    "path": path.as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )

    print(f"Wrote {MANIFEST_PATH} for {len(artefacts)} artefacts")


if __name__ == "__main__":
    main()
