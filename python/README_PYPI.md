# crispr

CUDA-accelerated CRISPR off-target search with MIT/CFD scoring. Includes a C++ core, CUDA kernels, CLI, and Python bindings.

Features:
- GPU/CPU parity for Hamming, MIT, CFD scoring (configurable via JSON tables).
- Build/load genome indexes; score guides with NGG/SpCas9 (20nt) v0.1.0 scope.
- CLI (`crispr-gpu`) and Python API (`import crispr_gpu as cg`).

Install:
```bash
pip install crispr
```

Source/build: see https://github.com/omniscoder/crispr-gpu
