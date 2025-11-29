# Benchmarks

This folder contains reproducible scripts to compare crispr-gpu to Cas-OFFinder and CRISPRitz on a small dataset (chr1, NGG, SpCas9, mismatches only).

## Data
- Use a single-chromosome FASTA (e.g., `hg38_chr1.fa`).
- Guide sets: `benchmarks/data/guides_chr1_10.tsv`, `guides_chr1_100.tsv`, `guides_chr1_1000.tsv` (provide your own or generate).

## Running
Example (crispr-gpu):
```
python benchmarks/scripts/run_crispr_gpu.py \
  --fasta benchmarks/data/hg38_chr1.fa \
  --guides benchmarks/data/guides_chr1_100.tsv \
  --out results_crispr_gpu.json
```

Cas-OFFinder (path to binary required):
```
python benchmarks/scripts/run_cas_offinder.py --cas-offinder /path/to/cas-offinder \
  --fasta benchmarks/data/hg38_chr1.fa \
  --guides benchmarks/data/guides_chr1_100.tsv \
  --out results_cas_offinder.json
```

CRISPRitz similarly:
```
python benchmarks/scripts/run_crispritz.py --crispritz /path/to/crispritz \
  --fasta benchmarks/data/hg38_chr1.fa \
  --guides benchmarks/data/guides_chr1_100.tsv \
  --out results_crispritz.json
```

Summarize:
```
python benchmarks/scripts/summarize.py --inputs results_*.json
```

## Report
Include hardware, CUDA version, and command lines. Suggested metrics: runtime, index time, hits/sec.
