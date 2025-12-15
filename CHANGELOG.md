# Changelog

## Compatibility
- Schema contracts are versioned and declared in outputs (`schema`, `schema_version`).
- As of `v0.2.0` (2025-12-14), the following schemas are frozen at v1: `report.v1`, `score_result.v1`, `bench_run.v1`, `hit.v1`.
- “Additive-only within v1” means: new optional fields may be added, but existing fields are never removed/renamed and their semantics do not change; breaking changes land as v2 side-by-side.

## 0.2.1 - 2025-12-15
- Supply chain:
  - Cosign keyless signing by digest for GHCR images
  - SPDX SBOM attestation published alongside the image
  - Publish workflow prints the digest and signature reference for verification
- Provenance:
  - `report.json` captures container digest and build metadata via `CRISPR_GPU_*` environment variables
  - Images bake build metadata so provenance is not best-effort
- Verification and citation rules:
  - `docs/supply_chain.md` shows how to verify signatures and SBOM attestations
  - Rule: do not cite `:latest`; cite a release tag or immutable digest
- Release surface:
  - Release assets auto-attach for future tags
  - README one-liner runs demo and validates schemas via `scripts/validate_artifacts.py`

## 0.2.0 - 2025-12-14
- Frozen, versioned output contracts: `report.v1`, `score_result.v1`, `bench_run.v1`, `hit.v1`
- Deterministic artifact bundle runner: `scripts/artifact_run.py` (demo + benchmarks + report)
- Docker images (CPU + CUDA) for artifact evaluation with predictable `/out/report.json`

## 0.1.0 - 2025-11-29
- CUDA and CPU parity engine with Hamming/MIT/CFD scoring
- Configurable scoring tables via JSON; default CFD/MIT tables bundled
- CLI and Python bindings with backend selection
- Base docs, README quickstart, and tests (C++ & Python)
