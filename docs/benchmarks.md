# Benchmarks (v0.1-benchmarked)

**Metric: CRISPR-GPU candidate throughput (CGCT)**  
Defined as: total candidate sites evaluated ÷ wall-clock time, for NGG PAM, guide length 20, K=4, Hamming score.

Hardware (local dev run):  
CPU: 12-core x86_64 (Ubuntu runner)  
GPU: NVIDIA GTX 1060 (6 GB, SM 6.1)  
Build: Release, CUDA on, warm GPU where noted.

## Synthetic genome (5 Mb, 624,487 sites), 50 guides
| Backend | Warmup | Time (s) | CGCT (candidates/s) | Hits |
| --- | --- | --- | --- | --- |
| CPU | n/a | 0.30 | ~1.04e8 | 8 |
| GPU | cold | 1.48 | ~2.11e7 | 8 |
| GPU | warm | 0.53 | ~5.9e7 | 8 |

## Synthetic genome (50 Mb), 50 guides
| Backend | Warmup | Time (s) | CGCT (candidates/s) | Hits |
| --- | --- | --- | --- | --- |
| CPU | n/a | 2.80 | ~1.12e8 | 114 |
| GPU | warm | 1.60 | ~1.95e8 | 114 |

Notes:
- Warm GPU runs use `CRISPR_GPU_WARMUP=1` to pay CUDA context cost before timing.
- Candidate count = (#sites) × (#guides); hits counted from output rows (header excluded).

## How to reproduce
```bash
cmake -B build -S . -DCRISPR_GPU_ENABLE_CUDA=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Small, cold:
./benchmarks/run_synthetic.sh

# Small, warm GPU:
CRISPR_GPU_WARMUP=1 ./benchmarks/run_synthetic.sh

# Large (50 Mb), warm GPU:
BENCH_SCALE=large CRISPR_GPU_WARMUP=1 ./benchmarks/run_synthetic.sh

# Guide sweep (50 and 500 guides):
GUIDE_SWEEP=50,500 ./benchmarks/run_synthetic.sh
```

## Kernel microbench (device-only)
Standalone device benchmark for the Hamming kernel (no host/index overhead):
```bash
cmake --build build -j --target kernel_microbench
MB_N=10000000 MB_ITERS=50 ./build/kernel_microbench
```
Outputs `kernel_candidates_per_sec` for direct GPU throughput inspection.
