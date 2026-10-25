#!/usr/bin/env python3
"""Run the complete single-cell RNA-seq immune profiling pipeline.

Usage:
    python run_pipeline.py          # Run all steps
    python run_pipeline.py --from 3 # Resume from step 3
"""

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = PROJECT_DIR / "scripts"

STEPS = [
    ("01 Load & QC", "01_load_and_qc.py"),
    ("02 Preprocess", "02_preprocess.py"),
    ("03 Reduce dimensions", "03_reduce_dimensions.py"),
    ("04 Cluster", "04_cluster.py"),
    ("05 Annotate cell types", "05_annotate_cell_types.py"),
    ("06 Trajectory inference", "06_trajectory.py"),
    ("07 T cell subclustering", "07_t_cell_subclustering.py"),
    ("08 Publication figures", "08_publication_figures.py"),
    ("09 Output manifest", "09_output_manifest.py"),
]

N_STEPS = len(STEPS)


def write_run_record(record):
    """Keep a machine-readable record, including failed and resumed runs."""
    results = PROJECT_DIR / "results"
    results.mkdir(exist_ok=True)
    (results / "run_metadata.json").write_text(json.dumps(record, indent=2) + "\n")


def git_output(*args):
    try:
        result = subprocess.run(
            ["git", *args], cwd=PROJECT_DIR, capture_output=True, text=True, check=False,
        )
    except OSError:
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def run_pipeline(start_from=1):
    if not 1 <= start_from <= N_STEPS:
        raise ValueError(f"Start step must be between 1 and {N_STEPS}.")

    print("=" * 60)
    print("Single-Cell RNA-seq Immune Profiling Pipeline")
    print("=" * 60)

    total_start = time.perf_counter()
    revision = git_output("rev-parse", "HEAD")
    changes = git_output("status", "--porcelain")
    record = {
        "started_at": datetime.now(timezone.utc).isoformat(),
        "status": "running", "start_step": start_from, "completed_steps": [],
        "python": platform.python_version(), "platform": platform.platform(),
        "git_commit": revision,
        "working_tree_dirty": bool(changes) if changes is not None else None,
        "packages": {dist.metadata["Name"]: dist.version for dist in importlib.metadata.distributions()},
    }
    write_run_record(record)
    for i, (name, script) in enumerate(STEPS, 1):
        if i < start_from:
            print(f"\n[{i}/{N_STEPS}] {name} -- SKIPPED")
            continue
        print(f"\n{'─' * 60}")
        print(f"[{i}/{N_STEPS}] {name}")
        print("─" * 60)
        step_start = time.perf_counter()
        result = subprocess.run(
            [sys.executable, str(SCRIPTS_DIR / script)],
            cwd=PROJECT_DIR,
            check=False,
        )
        if result.returncode != 0:
            record.update(status="failed", failed_step=i, exit_code=result.returncode)
            write_run_record(record)
            print(f"\nERROR: Step {i} failed with exit code {result.returncode}")
            sys.exit(result.returncode)
        elapsed = time.perf_counter() - step_start
        record["completed_steps"].append({"step": i, "script": script, "seconds": round(elapsed, 3)})
        write_run_record(record)
        print(f"  Completed in {elapsed:.1f}s")

    total = time.perf_counter() - total_start
    raw_input = PROJECT_DIR / "data" / "pbmc3k_raw.h5ad"
    if raw_input.exists():
        digest = hashlib.sha256()
        with raw_input.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        record["input_sha256"] = digest.hexdigest()
    record.update(status="complete", completed_at=datetime.now(timezone.utc).isoformat(),
                  seconds=round(total, 3))
    write_run_record(record)
    print(f"\n{'=' * 60}")
    print(f"Pipeline complete in {total:.1f}s")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Run single-cell RNA-seq pipeline")
    parser.add_argument("--from", dest="start_from", type=int, choices=range(1, N_STEPS + 1), default=1,
                        help=f"Step number to start from (1-{N_STEPS})")
    args = parser.parse_args()
    run_pipeline(start_from=args.start_from)


if __name__ == "__main__":
    main()
