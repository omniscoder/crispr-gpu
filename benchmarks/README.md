# Benchmarks

Minimal, reproducible GPU vs CPU/off-the-shelf comparisons. Scripts expect a single-chromosome FASTA (e.g., hg38 chr1) and small guide sets.

See docs/benchmarks.md for methodology.

Quick runs:
```bash
cmake -B build -S . -DCRISPR_GPU_ENABLE_CUDA=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./benchmarks/run_synthetic.sh
BENCH_JSONL=synthetic_runs.jsonl ./benchmarks/run_synthetic.sh
```

Kernel-only microbench:
```bash
cmake --build build -j --target kernel_microbench
./build/kernel_microbench --format kv
```
