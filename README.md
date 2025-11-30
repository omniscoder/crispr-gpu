# crispr-gpu

CUDA-accelerated CRISPR off-target search with configurable MIT/CFD scoring, C++/Python APIs, and a single CLI.

## Quickstart (CLI)

```bash
# Build index
crispr-gpu index \
  --fasta hg38.fa \
  --pam NGG \
  --guide-length 20 \
  --out hg38_spcas9_ngg.idx

# Score guides
crispr-gpu score \
  --index hg38_spcas9_ngg.idx \
  --guides guides.tsv \
  --max-mm 4 \
  --score-model cfd \
  --output hits.tsv
```

## Quickstart (Python)

```python
import crispr_gpu as cg

idx = cg.GenomeIndex.load("hg38_spcas9_ngg.idx")

params = cg.EngineParams()
params.score_params.model = cg.ScoreModel.CFD
params.max_mismatches = 4
params.backend = cg.Backend.GPU

engine = cg.OffTargetEngine(idx, params)

hits = engine.score_guide(
    cg.Guide(name="my_g1", sequence="GGGAAACCCGGGAAACCCGG", pam="NGG")
)
for h in hits[:5]:
    print(h.chrom_id, h.pos, h.mismatches, h.score)
```

## Install (from source)

```bash
pip install crispr            # from PyPI
# or build locally
cmake -B build -S . -DCRISPR_GPU_ENABLE_CUDA=ON -DCRISPR_GPU_BUILD_PYTHON=ON
cmake --build build --config Release
ctest --output-on-failure --test-dir build   # optional
pip install .
```

## Docs
- docs/getting_started.md
- docs/cfd_tables.md
- docs/benchmarks.md (skeleton)

## Synthetic Benchmark (quick sanity)
Synthetic genome, NGG, guide length 20, K=4, Hamming, 50 random guides.

| Backend | GPU warmup | Genome size | Time (s) |
| --- | --- | --- | --- |
| CPU | n/a | 5 Mb | ~0.30 |
| GPU | cold (includes CUDA init) | 5 Mb | ~1.48 |
| GPU | warm (CRISPR_GPU_WARMUP=1) | 5 Mb | ~0.53 |
| CPU | n/a | 50 Mb | ~2.80 |
| GPU | warm | 50 Mb | ~1.60 |

Run it yourself:
```bash
cmake -B build -S . -DCRISPR_GPU_ENABLE_CUDA=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./benchmarks/run_synthetic.sh              # CPU + GPU (if available)
CRISPR_GPU_WARMUP=1 ./benchmarks/run_synthetic.sh   # warm GPU timing
BENCH_SCALE=large ./benchmarks/run_synthetic.sh     # 50 Mb genome
```

## Version
0.1.0

## License
MIT
