# Getting Started

## Prerequisites
- CUDA toolkit (optional for GPU path; disable with `-DCRISPR_GPU_ENABLE_CUDA=OFF`)
- CMake >= 3.18
- C++14 compiler

## Build
```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Run CLI
```bash
./build/crispr-gpu index --fasta hg38.fa --pam NGG --guide-length 20 --out hg38.idx
./build/crispr-gpu score --index hg38.idx --guides guides.tsv --max-mm 4 --score-model hamming --output hits.tsv
```

## Python
```bash
pip install .
python -c "import crispr_gpu as cg; print(cg)"
```
