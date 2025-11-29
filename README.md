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

## Version
0.1.0

## License
MIT
