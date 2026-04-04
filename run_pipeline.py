#!/usr/bin/env python3
"""Run the complete single-cell RNA-seq immune profiling pipeline.

Usage:
    python run_pipeline.py          # Run all steps
    python run_pipeline.py --from 3 # Resume from step 3
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent / "scripts"

STEPS = [
    ("01 Load & QC", "01_load_and_qc.py"),
    ("02 Preprocess", "02_preprocess.py"),
    ("03 Reduce dimensions", "03_reduce_dimensions.py"),
    ("04 Cluster", "04_cluster.py"),
    ("05 Annotate cell types", "05_annotate_cell_types.py"),
    ("06 Trajectory inference", "06_trajectory.py"),
    ("07 T cell subclustering", "07_t_cell_subclustering.py"),
    ("08 Publication figures", "08_publication_figures.py"),
]

N_STEPS = len(STEPS)


def run_pipeline(start_from=1):
    print("=" * 60)
    print("Single-Cell RNA-seq Immune Profiling Pipeline")
    print("=" * 60)

    total_start = time.time()
    for i, (name, script) in enumerate(STEPS, 1):
        if i < start_from:
            print(f"\n[{i}/{N_STEPS}] {name} -- SKIPPED")
            continue
        print(f"\n{'─' * 60}")
        print(f"[{i}/{N_STEPS}] {name}")
        print("─" * 60)
        step_start = time.time()
        result = subprocess.run(
            [sys.executable, str(SCRIPTS_DIR / script)],
            cwd=str(Path(__file__).parent),
        )
        if result.returncode != 0:
            print(f"\nERROR: Step {i} failed with exit code {result.returncode}")
            sys.exit(result.returncode)
        elapsed = time.time() - step_start
        print(f"  Completed in {elapsed:.1f}s")

    total = time.time() - total_start
    print(f"\n{'=' * 60}")
    print(f"Pipeline complete in {total:.1f}s")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Run single-cell RNA-seq pipeline")
    parser.add_argument("--from", dest="start_from", type=int, default=1,
                        help=f"Step number to start from (1-{N_STEPS})")
    args = parser.parse_args()
    run_pipeline(start_from=args.start_from)


if __name__ == "__main__":
    main()
