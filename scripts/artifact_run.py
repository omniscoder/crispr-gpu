#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import platform
import hashlib
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


REPO_ROOT = Path(__file__).resolve().parents[1]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def run_cmd(
    cmd: list[str],
    *,
    cwd: Optional[Path] = None,
    env: Optional[dict[str, str]] = None,
    timeout_sec: Optional[float] = None,
) -> tuple[int, str, str, float]:
    start = time.perf_counter()
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout_sec,
    )
    dt = time.perf_counter() - start
    return p.returncode, p.stdout, p.stderr, dt


def maybe_cmd_output(cmd: list[str], *, cwd: Optional[Path] = None) -> Optional[str]:
    try:
        rc, out, err, _ = run_cmd(cmd, cwd=cwd, timeout_sec=10)
    except Exception:
        return None
    if rc != 0:
        return None
    s = (out or err).strip()
    return s or None


def find_exe(name: str, *, candidates: list[Path]) -> str:
    for p in candidates:
        if p.exists() and os.access(p, os.X_OK):
            return str(p)
    w = shutil.which(name)
    if w:
        return w
    raise FileNotFoundError(f"Unable to find executable: {name}")


def collect_git_provenance(repo_root: Path) -> dict[str, Any]:
    if not (repo_root / ".git").exists():
        return {"present": False}
    head = maybe_cmd_output(["git", "rev-parse", "HEAD"], cwd=repo_root)
    describe = maybe_cmd_output(["git", "describe", "--tags", "--always", "--dirty"], cwd=repo_root)
    status = maybe_cmd_output(["git", "status", "--porcelain=v1"], cwd=repo_root)
    dirty = bool(status and status.strip())
    return {
        "present": True,
        "head": head,
        "describe": describe,
        "dirty": dirty,
    }


def collect_system_provenance() -> dict[str, Any]:
    uname = platform.uname()
    return {
        "platform": platform.platform(),
        "uname": {
            "system": uname.system,
            "node": uname.node,
            "release": uname.release,
            "version": uname.version,
            "machine": uname.machine,
            "processor": uname.processor,
        },
        "python": {
            "version": sys.version.splitlines()[0],
            "executable": sys.executable,
        },
    }


def collect_tool_provenance(repo_root: Path) -> dict[str, Any]:
    tools: dict[str, Any] = {}
    tools["cmake"] = maybe_cmd_output(["cmake", "--version"])
    tools["cxx"] = maybe_cmd_output(["c++", "--version"]) or maybe_cmd_output(["g++", "--version"])
    tools["nvcc"] = maybe_cmd_output(["nvcc", "--version"])
    tools["nvidia_smi"] = maybe_cmd_output(["nvidia-smi"])
    # Relevant environment knobs (avoid dumping all env).
    tools["env"] = {k: os.environ[k] for k in sorted(os.environ.keys()) if k.startswith("CRISPR_GPU_")}
    return tools


def cuda_available_py() -> Optional[bool]:
    try:
        import crispr_gpu as cg  # type: ignore

        return bool(cg.cuda_available())
    except Exception:
        return None


@dataclass(frozen=True)
class DemoArtifacts:
    index_path: Path
    guides_path: Path
    cpu_json: Path
    gpu_json: Optional[Path]
    timings_sec: dict[str, float]


def run_demo(
    *,
    crispr_gpu: str,
    out_dir: Path,
    gpu_ok: bool,
) -> DemoArtifacts:
    demo_dir = out_dir / "demo"
    demo_dir.mkdir(parents=True, exist_ok=True)

    fasta = demo_dir / "demo.fa"
    guides = demo_dir / "demo_guides.tsv"
    index = demo_dir / "demo.idx"

    # Micro "realistic" demo: 20nt SpCas9-style guide with one engineered 1-mismatch off-target.
    # Use --plus-only to avoid reverse-strand PAM matches inflating the tiny example.
    guide = "GATCTACGATCTACGATCTA"  # 20nt, no "GG" runs (keeps PAM matches minimal)
    off1 = "GATCTACGATCTACGATCTC"   # 1 mismatch vs guide (last base A->C)
    pam1 = "AGG"  # NGG
    pam2 = "TGG"  # NGG
    spacer = "ATATATATATATATATATATATATATATATATATATATAT"
    seq = spacer + guide + pam1 + "A" + spacer + off1 + pam2 + "A" + spacer
    fasta.write_text(">chr1\n" + seq + "\n", encoding="utf-8")
    guides.write_text("name\tsequence\tpam\ng1\t" + guide + "\tNGG\n", encoding="utf-8")

    timings: dict[str, float] = {}

    rc, out, err, dt = run_cmd(
        [
            crispr_gpu,
            "index",
            "--fasta",
            str(fasta),
            "--pam",
            "NGG",
            "--guide-length",
            "20",
            "--plus-only",
            "--out",
            str(index),
        ],
        cwd=REPO_ROOT,
    )
    timings["index_build"] = dt
    (demo_dir / "index.stdout.txt").write_text(out, encoding="utf-8")
    (demo_dir / "index.stderr.txt").write_text(err, encoding="utf-8")
    if rc != 0:
        raise RuntimeError("demo index failed")

    cpu_json = demo_dir / "score_cpu.json"
    rc, out, err, dt = run_cmd(
        [
            crispr_gpu,
            "score",
            "--index",
            str(index),
            "--guides",
            str(guides),
            "--max-mm",
            "1",
            "--score-model",
            "hamming",
            "--backend",
            "cpu",
            "--search-backend",
            "brute",
            "--output-format",
            "json",
            "--sort",
            "--output",
            str(cpu_json),
        ],
        cwd=REPO_ROOT,
    )
    timings["score_cpu"] = dt
    (demo_dir / "score_cpu.stdout.txt").write_text(out, encoding="utf-8")
    (demo_dir / "score_cpu.stderr.txt").write_text(err, encoding="utf-8")
    if rc != 0:
        raise RuntimeError("demo score cpu failed")

    gpu_json: Optional[Path] = None
    if gpu_ok:
        gpu_json = demo_dir / "score_gpu.json"
        rc, out, err, dt = run_cmd(
            [
                crispr_gpu,
                "score",
                "--index",
                str(index),
                "--guides",
                str(guides),
                "--max-mm",
                "1",
                "--score-model",
                "hamming",
                "--backend",
                "gpu",
                "--search-backend",
                "brute",
                "--output-format",
                "json",
                "--sort",
                "--output",
                str(gpu_json),
            ],
            cwd=REPO_ROOT,
        )
        timings["score_gpu"] = dt
        (demo_dir / "score_gpu.stdout.txt").write_text(out, encoding="utf-8")
        (demo_dir / "score_gpu.stderr.txt").write_text(err, encoding="utf-8")
        if rc != 0:
            raise RuntimeError("demo score gpu failed")

    return DemoArtifacts(
        index_path=index,
        guides_path=guides,
        cpu_json=cpu_json,
        gpu_json=gpu_json,
        timings_sec=timings,
    )


def run_synthetic_bench(
    *,
    out_dir: Path,
    gpu_ok: bool,
    scale: str,
    quick: bool,
) -> Path:
    bench_jsonl = out_dir / "bench_synthetic.jsonl"
    bench_jsonl.unlink(missing_ok=True)

    env = dict(os.environ)
    env["BENCH_JSONL"] = str(bench_jsonl)
    env["BENCH_SCALE"] = scale
    env["CRISPR_GPU_WARMUP"] = "1" if gpu_ok else "0"
    env["SKIP_GPU"] = "0" if gpu_ok else "1"
    if quick:
        # Keep small enough for CI/docker but large enough to avoid /usr/bin/time rounding to 0.00s.
        env["GENOME_LEN"] = "100000"
        env["GUIDE_COUNT"] = "10"

    script = REPO_ROOT / "benchmarks" / "run_synthetic.sh"
    rc, out, err, _ = run_cmd([str(script)], cwd=REPO_ROOT, env=env)
    (out_dir / "bench_synthetic.stdout.txt").write_text(out, encoding="utf-8")
    (out_dir / "bench_synthetic.stderr.txt").write_text(err, encoding="utf-8")
    if rc != 0:
        raise RuntimeError("synthetic benchmark failed")
    return bench_jsonl


def maybe_run_kernel_microbench(
    *,
    out_dir: Path,
    gpu_ok: bool,
) -> Optional[dict[str, Any]]:
    if not gpu_ok:
        return None

    candidates = [
        REPO_ROOT / "build_cuda" / "kernel_microbench",
        REPO_ROOT / "build" / "kernel_microbench",
    ]
    try:
        exe = find_exe("kernel_microbench", candidates=candidates)
    except FileNotFoundError:
        return None

    out_path = out_dir / "kernel_microbench.json"
    rc, out, err, _ = run_cmd([exe, "--format", "json", "--output", str(out_path), "--iters", "10", "--warmup", "2"])
    (out_dir / "kernel_microbench.stdout.txt").write_text(out, encoding="utf-8")
    (out_dir / "kernel_microbench.stderr.txt").write_text(err, encoding="utf-8")
    if rc != 0 or not out_path.exists():
        return None
    return json.loads(out_path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def collect_artifact_hashes(out_dir: Path) -> dict[str, Any]:
    # Hash everything under out_dir except transient stdout/stderr capture.
    # (Those are kept for debugging but are not treated as stable artifacts.)
    hashes: dict[str, Any] = {}
    for p in sorted(out_dir.rglob("*")):
        if p.is_dir():
            continue
        rel = p.relative_to(out_dir).as_posix()
        if rel == "report.json":
            continue
        if rel.endswith(".stdout.txt") or rel.endswith(".stderr.txt"):
            continue
        hashes[rel] = {"sha256": sha256_file(p), "size_bytes": p.stat().st_size}
    return hashes


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


def write_report_json(report: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_report_md(report: dict[str, Any], path: Path) -> None:
    demo = report.get("demo", {})
    prov = report.get("provenance", {})
    git = prov.get("git", {})
    sysinfo = prov.get("system", {})
    tools = prov.get("tools", {})

    lines: list[str] = []
    lines.append(f"# crispr-gpu artifact report (v1)")
    lines.append("")
    lines.append(f"- Generated (UTC): `{report.get('generated_at_utc')}`")
    if git.get("present"):
        lines.append(f"- Git: `{git.get('describe')}`")
    lines.append(f"- Platform: `{sysinfo.get('platform')}`")
    lines.append(f"- Python: `{sysinfo.get('python', {}).get('version')}`")
    if tools.get("nvcc"):
        lines.append("- NVCC: present")
    if tools.get("nvidia_smi"):
        lines.append("- NVIDIA SMI: present")
    lines.append("")

    lines.append("## Demo")
    lines.append("")
    if demo.get("index_path"):
        lines.append(f"- Index: `{demo.get('index_path')}`")
        lines.append(f"- Guides: `{demo.get('guides_path')}`")
        lines.append(f"- CPU result: `{demo.get('cpu_score_result')}`")
        if demo.get("gpu_score_result"):
            lines.append(f"- GPU result: `{demo.get('gpu_score_result')}`")
        lines.append(f"- Timings (sec): `{demo.get('timings_sec')}`")
    else:
        lines.append("- Not run")
    lines.append("")

    lines.append("## Benchmarks")
    lines.append("")
    benches = report.get("benchmarks", [])
    lines.append(f"- Runs: `{len(benches)}`")
    for r in benches[:8]:
        p = r.get("params", {})
        t = r.get("timing_sec", {})
        lines.append(
            f"- synthetic scale={p.get('scale')} guides={p.get('guides')} K={p.get('k')} "
            f"cpu={t.get('cpu')} gpu_warm={t.get('gpu_warm')}"
        )
    if len(benches) > 8:
        lines.append(f"- (truncated; see JSON for full list)")
    lines.append("")

    km = report.get("kernel_microbench")
    lines.append("## Kernel microbench")
    lines.append("")
    if km:
        res = km.get("results", {})
        lines.append(f"- CGCT (cand/s): `{res.get('cgct_candidates_per_sec')}`")
        lines.append(f"- Mean sec: `{res.get('mean_sec')}`")
    else:
        lines.append("- Not available")
    lines.append("")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="Run demo + benchmarks and generate a stable report.")
    ap.add_argument("--out-dir", default=str(REPO_ROOT / "reports" / "latest"))
    ap.add_argument("--mode", choices=["all", "demo", "bench"], default="all")
    ap.add_argument("--bench-scale", choices=["small", "large"], default="small")
    ap.add_argument("--quick", action="store_true", help="Use tiny synthetic inputs for fast CI/demo runs.")
    ap.add_argument("--skip-gpu", action="store_true")
    ap.add_argument(
        "--redact",
        action="store_true",
        help="Remove host-specific paths/hostnames from the report (useful for publishing sample reports).",
    )
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    crispr_gpu = find_exe(
        "crispr-gpu",
        candidates=[
            REPO_ROOT / "build" / "crispr-gpu",
            REPO_ROOT / "build_cuda" / "crispr-gpu",
        ],
    )

    gpu_available = bool(cuda_available_py() or False)
    gpu_ok = gpu_available and not args.skip_gpu

    demo_artifacts: Optional[DemoArtifacts] = None
    bench_rows: list[Any] = []
    kernel_microbench: Optional[dict[str, Any]] = None

    cpu_score: Optional[dict[str, Any]] = None
    gpu_score: Optional[dict[str, Any]] = None

    if args.mode in ("all", "demo"):
        demo_artifacts = run_demo(crispr_gpu=crispr_gpu, out_dir=out_dir, gpu_ok=gpu_ok)
        cpu_score = load_json(demo_artifacts.cpu_json)
        gpu_score = load_json(demo_artifacts.gpu_json) if demo_artifacts.gpu_json else None

    if args.mode in ("all", "bench"):
        bench_jsonl = run_synthetic_bench(
            out_dir=out_dir, gpu_ok=gpu_ok, scale=args.bench_scale, quick=args.quick
        )
        bench_rows = load_jsonl(bench_jsonl)
        kernel_microbench = maybe_run_kernel_microbench(out_dir=out_dir, gpu_ok=gpu_ok)

    def rel_to_out(p: Optional[Path]) -> Optional[str]:
        if p is None:
            return None
        try:
            return str(p.relative_to(out_dir))
        except Exception:
            return str(p)

    report: dict[str, Any] = {
        "schema": "crispr-gpu/report/v1",
        "schema_version": 1,
        "generated_at_utc": utc_now_iso(),
        "provenance": {
            "git": collect_git_provenance(REPO_ROOT),
            "system": collect_system_provenance(),
            "tools": collect_tool_provenance(REPO_ROOT),
        },
        "demo": {},
        "benchmarks": bench_rows,
        "kernel_microbench": kernel_microbench,
    }

    if demo_artifacts:
        report["demo"] = {
            "index_path": rel_to_out(demo_artifacts.index_path),
            "guides_path": rel_to_out(demo_artifacts.guides_path),
            "cpu_score_result": rel_to_out(demo_artifacts.cpu_json),
            "gpu_score_result": rel_to_out(demo_artifacts.gpu_json),
            "timings_sec": demo_artifacts.timings_sec,
            "summary": {
                "cpu_num_hits": (cpu_score or {}).get("summary", {}).get("num_hits"),
                "gpu_num_hits": (gpu_score or {}).get("summary", {}).get("num_hits") if gpu_score else None,
            },
        }
    report["provenance"]["tools"]["crispr_gpu"] = {
        "path": crispr_gpu,
        "version": maybe_cmd_output([crispr_gpu, "--version"], cwd=REPO_ROOT),
    }

    schemas_used: dict[tuple[str, int], None] = {}
    schemas_used[("crispr-gpu/report/v1", 1)] = None
    for row in bench_rows:
        try:
            schemas_used[(str(row["schema"]), int(row["schema_version"]))] = None
        except Exception:
            pass
    if cpu_score and "schema" in cpu_score and "schema_version" in cpu_score:
        schemas_used[(str(cpu_score["schema"]), int(cpu_score["schema_version"]))] = None
    if kernel_microbench and "schema" in kernel_microbench and "schema_version" in kernel_microbench:
        schemas_used[(str(kernel_microbench["schema"]), int(kernel_microbench["schema_version"]))] = None

    report["schemas_used"] = [
        {"schema": k[0], "schema_version": k[1]} for k in sorted(schemas_used.keys())
    ]

    if args.redact:
        sysinfo = report.get("provenance", {}).get("system", {})
        if isinstance(sysinfo, dict):
            uname = sysinfo.get("uname")
            if isinstance(uname, dict):
                uname.pop("node", None)
            py = sysinfo.get("python")
            if isinstance(py, dict):
                py.pop("executable", None)
        tools = report.get("provenance", {}).get("tools", {})
        if isinstance(tools, dict):
            cg = tools.get("crispr_gpu")
            if isinstance(cg, dict) and "path" in cg:
                cg["path"] = "crispr-gpu"

    report_json = out_dir / "report.json"
    report_md = out_dir / "report.md"
    write_report_json(report, report_json)
    write_report_md(report, report_md)

    report["artifacts"] = collect_artifact_hashes(out_dir)
    # Re-write JSON now that we have hashes.
    write_report_json(report, report_json)

    print(f"Wrote: {report_json}")
    print(f"Wrote: {report_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
