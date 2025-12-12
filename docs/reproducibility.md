# Reproducibility Guide

This project aims to be publishable as a research artifact. This document is the “runbook” for rebuilding and reproducing the benchmark numbers.

## Build (CPU-only)

```bash
cmake -B build -S . \
  -DCRISPR_GPU_ENABLE_CUDA=OFF \
  -DCRISPR_GPU_BUILD_PYTHON=ON \
  -DCRISPR_GPU_ENABLE_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --output-on-failure --test-dir build
```

## Build (CUDA)

```bash
cmake -B build -S . \
  -DCRISPR_GPU_ENABLE_CUDA=ON \
  -DCRISPR_GPU_BUILD_PYTHON=ON \
  -DCRISPR_GPU_ENABLE_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --output-on-failure --test-dir build
```

Notes:
- Set `-DCMAKE_CUDA_ARCHITECTURES=...` if you want a tight binary for your GPU(s).
- Use `./build/crispr-gpu warmup` (or `CRISPR_GPU_WARMUP=1` in scripts) to avoid timing CUDA context creation.

## Capture environment metadata

For papers, record at least:

```bash
uname -a
cmake --version
c++ --version || g++ --version || clang++ --version
python3 --version

# If using CUDA:
nvidia-smi
nvcc --version
```

## Synthetic benchmark

The synthetic benchmark generates a deterministic random genome and guide set (seed=0) and reports end-to-end timing.

```bash
./benchmarks/run_synthetic.sh
CRISPR_GPU_WARMUP=1 ./benchmarks/run_synthetic.sh
BENCH_SCALE=large CRISPR_GPU_WARMUP=1 ./benchmarks/run_synthetic.sh
```

Machine-readable output:

```bash
BENCH_JSONL=synthetic_runs.jsonl BENCH_SCALE=large CRISPR_GPU_WARMUP=1 ./benchmarks/run_synthetic.sh
```

`synthetic_runs.jsonl` is JSONL (one JSON object per configuration).

## Kernel microbenchmark (device-only)

This isolates the candidate-scan kernel (single-guide scan) and reports `cgct_candidates_per_sec`.

```bash
cmake --build build -j --target kernel_microbench
./build/kernel_microbench --format kv
./build/kernel_microbench --format json --output kernel_microbench.json
```

Knobs:
- `--hit-fraction`: controls how many candidates pass the mismatch filter (affects atomic/write pressure).
- `--iters` / `--warmup`: stabilize timing for reporting.

