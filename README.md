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

# Score guides (Hamming, brute-force search)
crispr-gpu score \
  --index hg38_spcas9_ngg.idx \
  --guides guides.tsv \
  --max-mm 4 \
  --score-model cfd \
  --output hits.tsv

# Use FM-index backend (K=0 fast; K>0 experimental)
crispr-gpu score \
  --index hg38_spcas9_ngg.idx \
  --guides guides.tsv \
  --search-backend fmi \
  --max-mm 0 \
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
# Optional: FM backend (K=0 recommended publicly)
# params.search_backend = cg.SearchBackend.FMIndex

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
### Scoring modes

The scorer supports Hamming (default), MIT, and CFD.

CLI:

```bash
# Hamming (default)
crispr-gpu score --index hg38.idx --guides guides.tsv

# CFD scoring with bundled defaults
crispr-gpu score --index hg38.idx --guides guides.tsv --score-model cfd

# MIT scoring with a custom table
crispr-gpu score --index hg38.idx --guides guides.tsv \
  --score-model mit \
  --score-table data/mit_custom.json
```

Python:

```python
from crispr_gpu import OffTargetEngine, EngineParams, ScoreParams, SearchBackend, GenomeIndex

idx = GenomeIndex.load("hg38.idx")
params = EngineParams(
    search_backend=SearchBackend.FMIndex,
    score_params=ScoreParams(model='cfd', table_path='data/cfd_default.json'),
)
eng = OffTargetEngine(idx, params)
hits = eng.score_guides(guides)
```

Bundled defaults live under `data/cfd_default.json` and `data/mit_default.json`; supply `--score-table` / `table_path` to override.
