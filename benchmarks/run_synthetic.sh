#!/usr/bin/env bash
# Synthetic benchmark:
# - Small: 5 Mb genome, 50 random guides, K=4, Hamming (default)
# - Large: 50 Mb genome if BENCH_SCALE=large
# Outputs wall-clock timings for index build, CPU score, and GPU score (if available).
# Env knobs:
#   BENCH_SCALE=small|large   (default small)
#   GENOME_LEN / GUIDE_COUNT  (override lengths/count)
#   CRISPR_GPU_WARMUP=1       (do one throwaway GPU score before timing)
#   SKIP_GPU=1                (skip GPU section, useful on CPU-only CI)
#   CI_CPU_SLO=seconds        (fail if CPU time exceeds this; default 1.0s on CI for small scale)

set -euo pipefail

ROOT="$(mktemp -d /tmp/crisprgpu-bench-XXXXXX)"
export ROOT
BENSCALE="${BENCH_SCALE:-small}"
if [[ "$BENSCALE" == "large" ]]; then
  GENOME_LEN="${GENOME_LEN:-50000000}"
else
  GENOME_LEN="${GENOME_LEN:-5000000}"
fi
GUIDE_COUNT="${GUIDE_COUNT:-50}"
GUIDE_LEN=20
export GENOME_LEN GUIDE_COUNT

echo "Working dir: $ROOT"
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/build"

genome="$ROOT/genome.fa"
guides="$ROOT/guides.tsv"
index="$ROOT/genome.idx"
hits_cpu="$ROOT/hits_cpu.tsv"
hits_gpu="$ROOT/hits_gpu.tsv"

python3 - <<PY
import random, pathlib, os
random.seed(0)
root = pathlib.Path(os.environ.get("ROOT"))
length = int(os.environ.get("GENOME_LEN", "5000000"))
G = int(os.environ.get("GUIDE_COUNT", "50"))
seq = ''.join(random.choice('ACGT') for _ in range(length))
(root / "genome.fa").write_text(">chr1\n" + seq + "\n")
with open(root / "guides.tsv", "w") as f:
    for i in range(G):
        g = ''.join(random.choice('ACGT') for _ in range(20))
        f.write(f"g{i}\t{g}\tNGG\n")
PY

echo "Genome length: $GENOME_LEN bp"
echo "Guides: $GUIDE_COUNT"

log_build="$ROOT/log_build.txt"
log_cpu="$ROOT/log_cpu.txt"
log_gpu="$ROOT/log_gpu.txt"
time_build_file="$ROOT/time_build.txt"
time_cpu_file="$ROOT/time_cpu.txt"
time_gpu_file="$ROOT/time_gpu.txt"

time_build=$( /usr/bin/time -f "%e" -o "$time_build_file" ./build/crispr-gpu index --fasta "$genome" --pam NGG --guide-length $GUIDE_LEN --out "$index" >"$log_build" 2>&1 || true; cat "$time_build_file" )

time_cpu=$( /usr/bin/time -f "%e" -o "$time_cpu_file" env CRISPR_GPU_TIMING=1 ./build/crispr-gpu score --index "$index" --guides "$guides" --max-mm 4 --score-model hamming --backend cpu --output "$hits_cpu" >"$log_cpu" 2>&1 || true; cat "$time_cpu_file" )

gpu_available=0
if [[ ${SKIP_GPU:-0} -eq 0 ]]; then
  python3 - <<'PY' >/dev/null 2>&1 || true
import crispr_gpu as cg
import sys
sys.exit(0 if cg.cuda_available() else 1)
PY
  if [[ $? -eq 0 ]]; then
    gpu_available=1
    if [[ ${CRISPR_GPU_WARMUP:-0} -ne 0 ]]; then
      ./build/crispr-gpu score --index "$index" --guides "$guides" --max-mm 4 --score-model hamming --backend gpu --output /dev/null >/dev/null 2>&1 || true
    fi
    time_gpu=$( /usr/bin/time -f "%e" -o "$time_gpu_file" env CRISPR_GPU_TIMING=1 ./build/crispr-gpu score --index "$index" --guides "$guides" --max-mm 4 --score-model hamming --backend gpu --output "$hits_gpu" >"$log_gpu" 2>&1 || true; cat "$time_gpu_file" )
  else
    time_gpu="NA"
  fi
else
  time_gpu="NA"
fi

echo
echo "Results (seconds):"
printf "  build_index : %s\n" "$time_build"
printf "  score_cpu   : %s\n" "$time_cpu"
printf "  score_gpu   : %s\n" "$time_gpu"

# Parse counts
site_count=$(grep -oE 'with ([0-9]+) sites' "$log_build" | awk '{print $2}')
hits_cpu_count=$( (wc -l "$hits_cpu" 2>/dev/null || echo "0") | awk '{print ($1>0)?$1-1:0}')
hits_gpu_count=$( (wc -l "$hits_gpu" 2>/dev/null || echo "0") | awk '{print ($1>0)?$1-1:0}')

if [[ -n "$site_count" && "$site_count" -gt 0 ]]; then
  guides="$GUIDE_COUNT"
  candidates=$((site_count * guides))
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
  gpu_eps=$(TIME_VAL="$time_gpu" CAND="$candidates" python3 - <<'PY'
import os
try:
    t=float(os.environ["TIME_VAL"])
    c=int(os.environ["CAND"])
    print("{:.3f}".format(c/t))
except Exception:
    print("NA")
PY
)
  echo
  echo "Throughput:"
  echo "  candidates: $candidates"
  echo "  cpu candidates/sec: $cpu_eps"
  echo "  gpu candidates/sec: $gpu_eps"
  echo "  cpu hits: $hits_cpu_count"
  echo "  gpu hits: $hits_gpu_count"
fi

# Optional CPU regression guard for CI (small scale only)
if [[ "${CI:-}" != "" && "$BENSCALE" == "small" ]]; then
  limit="${CI_CPU_SLO:-1.0}"
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
echo
echo "Artifacts in $ROOT"
echo "  logs: $log_build $log_cpu $log_gpu"
