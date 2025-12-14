#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[Any]:
    rows: list[Any] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append(json.loads(line))
    return rows


def resolve_artifact_path(out_dir: Path, p: str) -> Path:
    cand = Path(p)
    if cand.is_absolute():
        return cand
    return out_dir / cand


def schema_file_for(schema_id: str, schema_version: int) -> Path:
    # Map "crispr-gpu/<name>/v1" -> "schemas/<name>.v1.json"
    prefix = "crispr-gpu/"
    if not schema_id.startswith(prefix) or "/v" not in schema_id:
        raise ValueError(f"Unsupported schema id: {schema_id}")
    name = schema_id[len(prefix) :].split("/v", 1)[0]
    return REPO_ROOT / "schemas" / f"{name}.v{schema_version}.json"


def require_jsonschema():
    try:
        import jsonschema  # type: ignore

        return jsonschema
    except Exception:
        print(
            "Missing dependency: jsonschema\n"
            "Install it with: python3 -m pip install jsonschema",
            file=sys.stderr,
        )
        raise SystemExit(2)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Validate an artifact_run.py output directory against bundled JSON Schemas."
    )
    ap.add_argument("out_dir", help="Directory containing report.json")
    args = ap.parse_args()

    jsonschema = require_jsonschema()

    out_dir = Path(args.out_dir)
    report_path = out_dir / "report.json"
    if not report_path.exists():
        raise SystemExit(f"Missing {report_path}")

    report = load_json(report_path)
    if report.get("schema") != "crispr-gpu/report/v1" or int(report.get("schema_version", 0)) != 1:
        raise SystemExit("Unexpected report schema header; expected crispr-gpu/report/v1 schema_version=1")

    report_schema = load_json(REPO_ROOT / "schemas" / "report.v1.json")
    jsonschema.validate(instance=report, schema=report_schema)

    # Validate demo artifacts (if present).
    demo = report.get("demo")
    if isinstance(demo, dict):
        cpu_p = demo.get("cpu_score_result")
        if isinstance(cpu_p, str):
            cpu_score_path = resolve_artifact_path(out_dir, cpu_p)
            cpu_score = load_json(cpu_score_path)
            score_schema = load_json(REPO_ROOT / "schemas" / "score_result.v1.json")
            jsonschema.validate(instance=cpu_score, schema=score_schema)

    # Validate benchmark rows (if present).
    bench_jsonl = out_dir / "bench_synthetic.jsonl"
    if bench_jsonl.exists():
        bench_rows = load_jsonl(bench_jsonl)
        bench_schema = load_json(REPO_ROOT / "schemas" / "bench_run.v1.json")
        for row in bench_rows:
            jsonschema.validate(instance=row, schema=bench_schema)

    # Validate kernel microbench (if present in report).
    km = report.get("kernel_microbench")
    if isinstance(km, dict) and km.get("schema") == "crispr-gpu/kernel_microbench/v1":
        km_schema = load_json(REPO_ROOT / "schemas" / "kernel_microbench.v1.json")
        jsonschema.validate(instance=km, schema=km_schema)

    # Ensure declared schemas resolve to bundled schema files (best-effort).
    schemas_used = report.get("schemas_used")
    if isinstance(schemas_used, list):
        for s in schemas_used:
            if not isinstance(s, dict):
                continue
            sid = s.get("schema")
            ver = s.get("schema_version")
            if not isinstance(sid, str) or not isinstance(ver, int):
                continue
            schema_path = schema_file_for(sid, ver)
            if not schema_path.exists():
                raise SystemExit(f"schemas_used references missing schema file: {schema_path}")

    print("ok: report.v1 (+ demo/bench if present) validated against bundled schemas")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

