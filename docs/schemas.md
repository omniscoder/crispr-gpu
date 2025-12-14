# JSON Schemas

This repo ships **versioned, stable JSON outputs** intended for tooling and artifact evaluation.

Schemas live in `schemas/`:

- `schemas/score_result.v1.json` — `crispr-gpu score --output-format json`
- `schemas/hit.v1.json` — `crispr-gpu score --output-format jsonl` (one hit per line)
- `schemas/bench_run.v1.json` — `benchmarks/run_synthetic.sh` with `BENCH_JSONL=...`
- `schemas/report.v1.json` — `scripts/artifact_run.py`
- `schemas/kernel_microbench.v1.json` — `kernel_microbench --format json`

Versioning policy:

- Schemas only change by adding `vN+1` (no breaking changes within a version).
- Each JSON includes `schema` and `schema_version`.
