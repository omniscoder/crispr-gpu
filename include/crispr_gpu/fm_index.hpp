#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <array>

namespace crispr_gpu {

// Minimal FM-index over DNA alphabet {A,C,G,T} with optional PAM filtering.
// Built over concatenated protospacers (guide_length) that satisfy PAM.
// This is a CPU-side structure for stage-1 candidate enumeration.

struct FmIndexHeader {
  uint8_t version{1};
  uint8_t guide_length{20};
  std::string pam{"NGG"};
  bool both_strands{true};
  uint64_t text_length{0};
  uint32_t occ_block{128};
  uint32_t sa_sample_rate{32};
};

struct FmIndex {
  FmIndexHeader header;
  // BWT over encoded alphabet {A=0,C=1,G=2,T=3,N=4}
  std::vector<uint8_t> bwt;
  // C array: cumulative counts for alphabet [A,C,G,T,N]. size=5.
  std::array<uint32_t,5> C{0,0,0,0,0};
  // Occ checkpoints: every occ_block rows, store counts for all symbols.
  std::vector<std::array<uint32_t,5>> occ_checkpoints;
  // Sampled SA values every sa_sample_rate rows.
  std::vector<uint32_t> sa_samples;
  // Text length (including sentinel); sentinel is the smallest symbol.
  uint32_t text_len{0};
  // Per-hit mapping back to contig coordinates.
  struct Locus {
    uint32_t chrom_id;
    uint32_t pos;
    uint8_t strand;
    uint32_t site_idx; // index into GenomeIndex sites_ for quick lookup
  };
  // One locus per suffix start; length = text_len.
  std::vector<Locus> loci;
};

// Build FM-index from a set of protospacers already filtered by PAM.
FmIndex build_fm_index(const std::vector<std::string> &protospacers,
                       const std::vector<FmIndex::Locus> &loci,
                       const FmIndexHeader &header);

// Serialize/deserialize.
void save_fm_index(const FmIndex &index, const std::string &path);
FmIndex load_fm_index(const std::string &path);

// Search allowing up to K mismatches (Hamming) via bounded DFS/backward search.
// Returns indices into `loci`.
std::vector<uint32_t> fm_search_hamming(const FmIndex &index,
                                        const std::string &pattern,
                                        uint8_t max_mismatches);

// Convert BWT row to reference position.
uint32_t fm_row_to_pos(const FmIndex &index, uint32_t row);


} // namespace crispr_gpu
