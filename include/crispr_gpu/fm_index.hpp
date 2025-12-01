#pragma once

#include <cstdint>
#include <string>
#include <vector>

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
  uint64_t sampled_sa_interval{32};
};

struct FmIndex {
  FmIndexHeader header;
  // C array: cumulative counts for alphabet [A,C,G,T]. size=5 for convenience.
  uint64_t C[5]{};
  // Occ counts sampled every occ_block for each symbol; flattened as [symbol][block].
  uint32_t occ_block{128};
  std::vector<uint32_t> occ_A;
  std::vector<uint32_t> occ_C;
  std::vector<uint32_t> occ_G;
  std::vector<uint32_t> occ_T;
  // Sampled suffix array (positions into concatenated text).
  std::vector<uint32_t> sa_samples;
  uint32_t sa_sample_rate{32};
  // Per-protospacer mapping back to chromosome/pos/strand.
  struct Locus {
    uint32_t chrom_id;
    uint32_t pos;
    uint8_t strand;
  };
  std::vector<Locus> loci; // same length as text_length/guide_length
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

} // namespace crispr_gpu

