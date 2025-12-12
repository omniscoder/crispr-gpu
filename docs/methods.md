# Methods / Design Notes

This document captures implementation-level details for reproducibility and review. It is intentionally “close to the code” (v0.1.0).

## Index construction

`GenomeIndex::build(fasta, params)` enumerates candidate protospacer sites from each FASTA contig:

- **PAM matching** uses IUPAC rules (`N`, `R`, `Y`, …). See `src/genome_index_build.cpp`.
- **Plus strand**: for each offset `i`, check `seq[i+guide_len : i+guide_len+pam_len]` against `pam`. If it matches, take the protospacer `seq[i : i+guide_len]`.
- **Minus strand (optional)**: check `seq[i : i+pam_len]` against `revcomp(pam)`. If it matches, take `seq[i+pam_len : i+pam_len+guide_len]`, reverse-complement it, and store as the protospacer sequence.
- **Ambiguous bases inside the protospacer** (anything other than A/C/G/T) are skipped because 2‑bit encoding rejects them.

Each accepted site is stored as a `SiteRecord`:

- `seq_bits`: 2‑bit packed protospacer sequence (guide orientation)
- `chrom_id`: contig index within the FASTA
- `pos`: 0‑based coordinate of the protospacer start on the forward reference sequence
- `strand`: `0` for `+`, `1` for `-`

## 2-bit encoding

`encode_sequence_2bit()` packs up to 32 bases into a `uint64_t` using `A=0, C=1, G=2, T=3`. The first base in the string becomes the most significant 2‑bit pair within the used region.

Hamming distance is computed using bitwise operations plus `popcount` over a length-derived mask (`hamming_distance_2bit()` in `src/types.cpp`).

## Candidate search backends

Stage‑1 candidate enumeration happens via `EngineParams.search_backend`:

- `brute` (production for mismatches): scan all sites in the index.
- `fmi` (experimental): per‑contig FM-index built over protospacers. The implementation supports Hamming‑K search internally, but the engine currently restricts this backend to **K=0 exact search** (see `src/engine.cpp`) to keep behavior conservative.

## Scoring models

Scoring is applied after mismatch filtering:

- `hamming`: simple monotone score `1/(1+mm)` (baseline sanity score).
- `mit` / `cfd`: lightweight per‑position/per‑type weights. Defaults are placeholders; they can be overridden via a JSON file (see `docs/cfd_tables.md`).

On GPU builds, scoring tables are uploaded once per process to device constant memory (`src/kernels.cu`).

## GPU pipeline (brute backend)

For GPU runs the engine:

1. Uploads `SiteRecord[]` to device memory (cached across calls).
2. Encodes the guide to 2‑bit form on the host.
3. Launches a single kernel (`src/kernels.cu`) that:
   - computes mismatches,
   - filters by `max_mismatches`,
   - computes score (Hamming/MIT/CFD),
   - appends passing hits via `atomicAdd` into a preallocated hit buffer.

The default kernel launch uses `block=256` and caps the 1‑D grid at `65535` blocks.

## Known limitations (v0.1.0)

- Index building currently reads whole contigs into memory (not streamed).
- FM backend is exact-only at the engine level (K=0) even though the FM search routine supports Hamming‑K.
- Default MIT/CFD tables are placeholders; provide validated weights for publications.

