# Changelog

## 0.2.0 - 2025-12-14
- Frozen, versioned output contracts: `report.v1`, `score_result.v1`, `bench_run.v1`, `hit.v1`
- Deterministic artifact bundle runner: `scripts/artifact_run.py` (demo + benchmarks + report)
- Docker images (CPU + CUDA) for artifact evaluation with predictable `/out/report.json`

## 0.1.0 - 2025-11-29
- CUDA and CPU parity engine with Hamming/MIT/CFD scoring
- Configurable scoring tables via JSON; default CFD/MIT tables bundled
- CLI and Python bindings with backend selection
- Base docs, README quickstart, and tests (C++ & Python)
