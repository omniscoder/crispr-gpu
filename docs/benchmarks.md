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
| CPU | n/a | 0.26 | 1.20e8 | 8 |
| GPU | cold | 0.86 | 3.63e7 | 8 |
| GPU | warm | 0.53 | 5.89e7 | 8 |

## Synthetic genome (50 Mb, 6,246,000 sites), 500 guides
| Backend | Warmup | Time (s) | CGCT (candidates/s) | Hits | Notes |
| --- | --- | --- | --- | --- | --- |
| CPU | n/a | 21.32 | 1.46e8 | 1281 | brute, Hamming |
| GPU | cold | 2.71 | 1.15e9 | 1281 | brute, Hamming |
| GPU | warm | 2.15 | 1.45e9 | 1281 | brute, Hamming |
| GPU | warm | 1.64 | 1.90e9 | 1282 | brute, Hamming, direct run |
| GPU | warm | 1.65 | 1.89e9 | 1282 | brute, CFD, direct run |

The last two rows come from direct GPU runs (no CPU pass) on the same index/guide set; they show the “brag” number for K=4 under Hamming and CFD scoring on this GTX 1060 (~1.9×10⁹ candidates/s).

### Guide sweep (50 Mb, warm GPU)
Candidate count scales as 6,246,000 sites × guides.

| Guides | gpu.stage1 (s) | gpu.stage2 (s) | gpu.guides (s) | gpu.guides CGCT |
| --- | --- | --- | --- | --- |
| 1   | ~0.0001 | ~0.58 | 1.35 | 4.63e6 |
| 10  | ~0.0001 | ~0.58 | 1.44 | 4.34e7 |
| 50  | ~0.0001 | ~0.58 | 1.40 | 2.23e8 |
| 500 | ~0.0001 | ~0.58 | 2.15 | 1.45e9 |

Notes:
- Warm GPU runs use `CRISPR_GPU_WARMUP=1` to pay CUDA context cost before timing.
- Candidate count = (#sites) × (#guides); hits counted from output rows (header excluded).
- Warm GPU 500-guide runs show ~±10% variance (1.29e9–1.55e9 cand/s) across repeated runs on the same GTX 1060; numbers above use a representative run.
- Persistent GPU state + batch scoring improved the 50 Mb / 500 guide warm CGCT from ~1.20e9 (v0.10-bench) to ~1.45e9 (~22% uplift) without kernel changes; this is the brute-force Hamming baseline for future FM-index and DPX work.
- FM backend: `--search-backend fmi` is supported. FM K=0 is modestly faster on synthetic large cases. **FM K>0 is disabled in public builds** (the engine throws if K>0 with FM); use the brute backend for mismatches.

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
./build/kernel_microbench
```
Configuration: grid 24,399, block 256, candidates 6.246e6, max_mm=4.
Latest result (GTX 1060): 2.01e9 candidates/s in 0.003107 s (single-guide scan of the 50 Mb index size).
Outputs `cgct_candidates_per_sec` for direct GPU throughput inspection.
