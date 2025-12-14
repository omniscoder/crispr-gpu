import copy
import json
import os
import subprocess
from pathlib import Path

import pytest


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path):
    return json.loads(path.read_text())


def _normalize_report_for_golden(doc: dict) -> dict:
    d = copy.deepcopy(doc)

    # Unstable-by-design fields.
    d.pop("generated_at_utc", None)

    # Machine-specific and commit-specific.
    if isinstance(d.get("provenance"), dict):
        d["provenance"] = {}

    # Per-run hashes will vary with per-run timings; ensure presence elsewhere.
    if "artifacts" in d:
        d["artifacts"] = {}

    # Demo timings are inherently machine-dependent.
    if isinstance(d.get("demo"), dict):
        d["demo"].pop("timings_sec", None)

    # Bench timings/throughput are machine-dependent.
    benches = d.get("benchmarks")
    if isinstance(benches, list):
        for r in benches:
            if not isinstance(r, dict):
                continue
            r.pop("generated_at_utc", None)
            if "timing_sec" in r:
                r["timing_sec"] = {}
            if "cgct_candidates_per_sec" in r:
                r["cgct_candidates_per_sec"] = {}

    return d


def test_artifact_run_golden_report(tmp_path):
    jsonschema = pytest.importorskip("jsonschema")

    repo_root = _repo_root()
    out_dir = tmp_path / "artifact"
    out_dir.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    env.update({"LC_ALL": "C", "LANG": "C", "PYTHONHASHSEED": "0"})

    subprocess.check_call(
        [
            "python3",
            "scripts/artifact_run.py",
            "--out-dir",
            str(out_dir),
            "--quick",
            "--skip-gpu",
            "--redact",
        ],
        cwd=repo_root,
        env=env,
    )

    report_path = out_dir / "report.json"
    assert report_path.exists()
    report = _load_json(report_path)

    # Basic sanity + schema headers.
    assert report["schema"] == "crispr-gpu/report/v1"
    assert report["schema_version"] == 1
    assert "schemas_used" in report

    # Validate key produced artifacts exist.
    demo = report.get("demo", {})
    assert isinstance(demo, dict)
    assert (out_dir / demo["cpu_score_result"]).exists()
    assert (out_dir / "bench_synthetic.jsonl").exists()

    # JSON Schema validation (best-effort).
    report_schema = _load_json(repo_root / "schemas" / "report.v1.json")
    bench_schema = _load_json(repo_root / "schemas" / "bench_run.v1.json")
    score_schema = _load_json(repo_root / "schemas" / "score_result.v1.json")
    jsonschema.validate(instance=report, schema=report_schema)

    cpu_score = _load_json(out_dir / demo["cpu_score_result"])
    jsonschema.validate(instance=cpu_score, schema=score_schema)

    bench_rows = []
    for line in (out_dir / "bench_synthetic.jsonl").read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        bench_rows.append(json.loads(line))
    assert bench_rows
    for row in bench_rows:
        jsonschema.validate(instance=row, schema=bench_schema)

    # Golden snapshot compare (normalized).
    golden_path = repo_root / "reports" / "sample" / "report.json"
    golden = _load_json(golden_path)

    assert _normalize_report_for_golden(report) == _normalize_report_for_golden(golden)

