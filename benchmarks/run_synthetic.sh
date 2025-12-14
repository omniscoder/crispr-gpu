#!/usr/bin/env bash
# Synthetic benchmark:
# - Small: 5 Mb genome, 50 random guides, K=4, Hamming (default)
# - Large: 50 Mb genome if BENCH_SCALE=large
# Outputs wall-clock timings for index build, CPU score, and GPU score (if available).
# Env knobs:
#   BENCH_SCALE=small|large   (default small)
#   GENOME_LEN / GUIDE_COUNT  (override lengths/count)
#   GUIDE_SWEEP=50,500,5000   (comma-separated counts; overrides GUIDE_COUNT loop)
#   CRISPR_GPU_WARMUP=1       (do one warm GPU pass after cold measurement)
#   SKIP_GPU=1                (skip GPU section, useful on CPU-only CI)
#   CI_CPU_SLO=seconds        (fail if CPU time exceeds this; default 1.0s on CI for small scale)
#   SEARCH_BACKEND=brute|fmi  (search backend; default brute)
#   K_SWEEP=0,1,2,4           (optional comma list of K values; overrides MAX_MM)
#   SCORE_MODEL=hamming|mit|cfd (default hamming)
#   BENCH_JSONL=path          (append one JSON object per run)

set -euo pipefail
export LC_ALL=C
export LANG=C
export PYTHONHASHSEED=0

ROOT="$(mktemp -d /tmp/crisprgpu-bench-XXXXXX)"
export ROOT
BENSCALE="${BENCH_SCALE:-small}"
if [[ "$BENSCALE" == "large" ]]; then
  GENOME_LEN="${GENOME_LEN:-50000000}"
else
  GENOME_LEN="${GENOME_LEN:-5000000}"
fi
GUIDE_COUNT="${GUIDE_COUNT:-50}"
if [[ -n "${GUIDE_SWEEP:-}" ]]; then
  IFS=',' read -ra GUIDE_LIST <<< "$GUIDE_SWEEP"
else
  GUIDE_LIST=("$GUIDE_COUNT")
fi
MAX_MM_VAL=${MAX_MM:-4}
# K sweep
if [[ -n "${K_SWEEP:-}" ]]; then
  IFS=',' read -ra K_LIST <<< "$K_SWEEP"
else
  K_LIST=(${MAX_MM_VAL})
fi
BACKEND="${SEARCH_BACKEND:-brute}"
SCORE_MODEL="${SCORE_MODEL:-hamming}"
# determine max guides needed
MAX_GUIDES="$GUIDE_COUNT"
for g in "${GUIDE_LIST[@]}"; do
  if [[ "$g" -gt "$MAX_GUIDES" ]]; then MAX_GUIDES="$g"; fi
done
export MAX_GUIDES
GUIDE_LEN=20
export GENOME_LEN GUIDE_COUNT

echo "Working dir: $ROOT"
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/build:$(pwd)/build/python"

CRISPR_GPU_CLI="${CRISPR_GPU_BIN:-}"
if [[ -z "$CRISPR_GPU_CLI" ]]; then
  if [[ -x "./build/crispr-gpu" ]]; then
    CRISPR_GPU_CLI="./build/crispr-gpu"
  else
    CRISPR_GPU_CLI="crispr-gpu"
  fi
fi

genome="$ROOT/genome.fa"
guides_max="$ROOT/guides_max.tsv"
index="$ROOT/genome.idx"

python3 - <<PY
import random, pathlib, os
SEED = 0
random.seed(SEED)
root = pathlib.Path(os.environ.get("ROOT"))
length = int(os.environ.get("GENOME_LEN", "5000000"))
G = int(os.environ.get("MAX_GUIDES", "50"))
seq = ''.join(random.choice('ACGT') for _ in range(length))
(root / "genome.fa").write_text(">chr1\n" + seq + "\n")
with open(root / "guides_max.tsv", "w") as f:
    for i in range(G):
        g = ''.join(random.choice('ACGT') for _ in range(20))
        f.write(f"g{i}\t{g}\tNGG\n")
PY

echo "Genome length: $GENOME_LEN bp"
echo "Guides: $GUIDE_COUNT"

log_build="$ROOT/log_build.txt"
time_build_file="$ROOT/time_build.txt"

time_build=$( /usr/bin/time -f "%e" -o "$time_build_file" "$CRISPR_GPU_CLI" index --fasta "$genome" --pam NGG --guide-length $GUIDE_LEN --out "$index" >"$log_build" 2>&1 || true; cat "$time_build_file" )

site_count=$(grep -oE 'with ([0-9]+) sites' "$log_build" | awk '{print $2}' || true)
site_count=${site_count:-0}

gpu_available=0
if [[ ${SKIP_GPU:-0} -eq 0 ]]; then
  python3 - <<'PY' >/dev/null 2>&1 || true
import crispr_gpu as cg
import sys
sys.exit(0 if cg.cuda_available() else 1)
PY
  [[ $? -eq 0 ]] && gpu_available=1
fi

echo
echo "Results (seconds and CGCT):"
printf "  build_index : %s\n" "$time_build"
echo

cold_gpu_logged=0

for guides_cur in "${GUIDE_LIST[@]}"; do
  guides_file="$ROOT/guides_${guides_cur}.tsv"
  head -n "$guides_cur" "$guides_max" > "$guides_file"

  for K in "${K_LIST[@]}"; do
    log_cpu="$ROOT/log_cpu_${guides_cur}_k${K}.txt"
    log_gpu="$ROOT/log_gpu_${guides_cur}_k${K}.txt"
    time_cpu_file="$ROOT/time_cpu_${guides_cur}_k${K}.txt"
    time_gpu_file="$ROOT/time_gpu_${guides_cur}_k${K}.txt"
    hits_cpu="$ROOT/hits_cpu_${guides_cur}_k${K}.tsv"
    hits_gpu="$ROOT/hits_gpu_${guides_cur}_k${K}.tsv"

    time_cpu=$( /usr/bin/time -f "%e" -o "$time_cpu_file" env CRISPR_GPU_TIMING=1 "$CRISPR_GPU_CLI" score --index "$index" --guides "$guides_file" --max-mm "$K" --score-model "$SCORE_MODEL" --backend cpu --search-backend "$BACKEND" --output "$hits_cpu" >"$log_cpu" 2>&1 || true; cat "$time_cpu_file" )
    time_gpu_cold="NA"
    time_gpu_warm="NA"

    if [[ $gpu_available -eq 1 ]]; then
      # cold run
      set +e
      /usr/bin/time -f "%e" -o "$time_gpu_file" env CRISPR_GPU_TIMING=1 "$CRISPR_GPU_CLI" score --index "$index" --guides "$guides_file" --max-mm "$K" --score-model "$SCORE_MODEL" --backend gpu --search-backend "$BACKEND" --output "$hits_gpu" >"$log_gpu" 2>&1
      rc=$?
      set -e
      if [[ $rc -eq 0 ]]; then
        time_gpu_cold=$(cat "$time_gpu_file")
      else
        time_gpu_cold="ERR"
      fi
      # warm run (optional)
      if [[ ${CRISPR_GPU_WARMUP:-0} -ne 0 ]]; then
        "$CRISPR_GPU_CLI" warmup >/dev/null 2>&1 || true
        set +e
        /usr/bin/time -f "%e" -o "$time_gpu_file" env CRISPR_GPU_TIMING=1 "$CRISPR_GPU_CLI" score --index "$index" --guides "$guides_file" --max-mm "$K" --score-model "$SCORE_MODEL" --backend gpu --search-backend "$BACKEND" --output "$hits_gpu" >"$log_gpu" 2>&1
        rc=$?
        set -e
        if [[ $rc -eq 0 ]]; then
          time_gpu_warm=$(cat "$time_gpu_file")
        else
          time_gpu_warm="ERR"
        fi
      else
        time_gpu_warm="NA"
      fi
    fi

    # stats
    hits_cpu_count=$( (wc -l "$hits_cpu" 2>/dev/null || echo "0") | awk '{print ($1>0)?$1-1:0}')
    hits_gpu_count=$( (wc -l "$hits_gpu" 2>/dev/null || echo "0") | awk '{print ($1>0)?$1-1:0}')
    candidates=$((site_count * guides_cur))

    cpu_eps=$(TIME_VAL="$time_cpu" CAND="$candidates" python3 - <<'PY'
import os
try:
    t=float(os.environ["TIME_VAL"])
    c=int(os.environ["CAND"])
    print("{:.3f}".format(c/t))
except Exception:
    print("NA")
PY
)
    gpu_eps_cold="NA"
    if [[ "$time_gpu_cold" != "NA" ]]; then
      gpu_eps_cold=$(TIME_VAL="$time_gpu_cold" CAND="$candidates" python3 - <<'PY'
import os
try:
    t=float(os.environ["TIME_VAL"])
    c=int(os.environ["CAND"])
    print("{:.3f}".format(c/t))
except Exception:
    print("NA")
PY
)
    fi

    gpu_eps_warm="NA"
    if [[ "$time_gpu_warm" != "NA" ]]; then
      gpu_eps_warm=$(TIME_VAL="$time_gpu_warm" CAND="$candidates" python3 - <<'PY'
import os
try:
    t=float(os.environ["TIME_VAL"])
    c=int(os.environ["CAND"])
    print("{:.3f}".format(c/t))
except Exception:
    print("NA")
PY
)
    fi

    echo "scale=$BENSCALE backend=$BACKEND score_model=$SCORE_MODEL K=$K guides=$guides_cur candidates=$candidates cpu_time=$time_cpu cpu_cgct=$cpu_eps gpu_cold_time=$time_gpu_cold gpu_cold_cgct=$gpu_eps_cold gpu_warm_time=$time_gpu_warm gpu_warm_cgct=$gpu_eps_warm hits_cpu=$hits_cpu_count hits_gpu=$hits_gpu_count"
    echo

    if [[ -n "${BENCH_JSONL:-}" ]]; then
      BENCH_JSONL_PATH="$BENCH_JSONL" \
      SCALE="$BENSCALE" SEARCH_BACKEND="$BACKEND" SCORE_MODEL="$SCORE_MODEL" \
      K="$K" GUIDES="$guides_cur" CANDIDATES="$candidates" SITE_COUNT="$site_count" \
      CPU_TIME="$time_cpu" CPU_CGCT="$cpu_eps" \
      GPU_COLD_TIME="$time_gpu_cold" GPU_COLD_CGCT="$gpu_eps_cold" \
      GPU_WARM_TIME="$time_gpu_warm" GPU_WARM_CGCT="$gpu_eps_warm" \
      HITS_CPU="$hits_cpu_count" HITS_GPU="$hits_gpu_count" \
      BUILD_INDEX_TIME="$time_build" \
      GENOME_LEN_BP="$GENOME_LEN" GUIDE_LEN="$GUIDE_LEN" SEED="0" \
      python3 - <<'PY'
import json, os
from datetime import datetime, timezone

def f(key):
    v = os.environ.get(key, "NA")
    if v in ("NA", "ERR", "", "None"):
        return None
    try:
        return float(v)
    except Exception:
        return None

row = {
    "schema": "crispr-gpu/bench_run/v1",
    "schema_version": 1,
    "tool": "crispr-gpu",
    "bench": "synthetic",
    "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
    "params": {
        "scale": os.environ.get("SCALE"),
        "genome_len_bp": int(os.environ.get("GENOME_LEN_BP", "0")),
        "guide_len": int(os.environ.get("GUIDE_LEN", "0")),
        "seed": int(os.environ.get("SEED", "0")),
        "guides": int(os.environ.get("GUIDES", "0")),
        "k": int(os.environ.get("K", "0")),
        "search_backend": os.environ.get("SEARCH_BACKEND"),
        "score_model": os.environ.get("SCORE_MODEL"),
    },
    "derived": {
        "site_count": int(os.environ.get("SITE_COUNT", "0")),
        "candidates": int(os.environ.get("CANDIDATES", "0")),
    },
    "timing_sec": {
        "build_index": f("BUILD_INDEX_TIME"),
        "cpu": f("CPU_TIME"),
        "gpu_cold": f("GPU_COLD_TIME"),
        "gpu_warm": f("GPU_WARM_TIME"),
    },
    "cgct_candidates_per_sec": {
        "cpu": f("CPU_CGCT"),
        "gpu_cold": f("GPU_COLD_CGCT"),
        "gpu_warm": f("GPU_WARM_CGCT"),
    },
    "hits": {
        "cpu": int(os.environ.get("HITS_CPU", "0")),
        "gpu": int(os.environ.get("HITS_GPU", "0")),
    },
}
with open(os.environ["BENCH_JSONL_PATH"], "a") as out:
    out.write(json.dumps(row, sort_keys=True) + "\n")
PY
    fi
  done
done

# Optional CPU regression guard for CI (small scale only, first CPU run)
# If SKIP_GPU=1 (CPU-only), either set CI_CPU_SLO explicitly or default to skipping the guard.
if [[ "${CI:-}" != "" && "$BENSCALE" == "small" ]]; then
  limit_env="${CI_CPU_SLO:-}"
  if [[ "$SKIP_GPU" == "1" && -z "$limit_env" ]]; then
    echo "[bench] SKIP_GPU=1 on CI with no CI_CPU_SLO set; skipping CPU SLO guard." >&2
  else
    limit="${limit_env:-1.0}"
    TIME_VAL="$time_cpu" LIMIT_VAL="$limit" python3 - <<'PY'
import os, sys
try:
    t=float(os.environ["TIME_VAL"])
    limit=float(os.environ["LIMIT_VAL"])
    if t > limit:
        print(f"CPU time {t:.3f}s exceeds limit {limit:.3f}s", file=sys.stderr)
        sys.exit(1)
except Exception:
    sys.exit(0)
PY
  fi
fi

echo "Artifacts in $ROOT"
echo "  logs: $log_build and per-guides logs under $ROOT"
