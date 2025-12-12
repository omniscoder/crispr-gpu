# Contributing

The goal is a research-grade, reproducible CRISPR off-target engine with GPU acceleration. Contributions that improve correctness, determinism, and measurement quality are especially welcome.

## Build and test

CPU-only:
```bash
cmake -B build -S . -DCRISPR_GPU_ENABLE_CUDA=OFF -DCRISPR_GPU_BUILD_PYTHON=ON -DCRISPR_GPU_ENABLE_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --output-on-failure --test-dir build
```

CUDA:
```bash
cmake -B build -S . -DCRISPR_GPU_ENABLE_CUDA=ON -DCRISPR_GPU_BUILD_PYTHON=ON -DCRISPR_GPU_ENABLE_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --output-on-failure --test-dir build
```

## Benchmarks

See `docs/benchmarks.md` and `docs/reproducibility.md`.

Quick sanity:
```bash
./benchmarks/run_synthetic.sh
BENCH_JSONL=synthetic_runs.jsonl ./benchmarks/run_synthetic.sh
```

Kernel-only:
```bash
cmake --build build -j --target kernel_microbench
./build/kernel_microbench --format kv
```

## What to include in PRs

- A clear motivation (bugfix, speedup, new scoring model, new search backend).
- Tests for new logic (C++ Catch2 and/or Python pytest).
- If performance-related: benchmark output and hardware/software metadata.

