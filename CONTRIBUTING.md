# Contributing

Contributions are welcome when they improve reproducibility, biological annotation, testing, or figure quality.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest -v
python run_pipeline.py
python scripts/validate_outputs.py
```

## Pull request checklist

- `pytest -v` passes.
- `ruff check` passes for changed Python files.
- `python scripts/validate_outputs.py` passes after any pipeline-output change.
- Changed figures are regenerated intentionally.
- Annotation changes document which marker genes or scoring rules changed.

## Scientific correctness

Cluster labels and trajectory claims should be evidence-backed. If you change cell-type markers, Leiden resolution selection, QC thresholds, or pseudotime rooting, include the effect on cell counts and key figures.

## Licence

By contributing, you agree that your contributions are licensed under the MIT licence.
