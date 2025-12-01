#pragma once

#include "crispr_gpu/types.hpp"
#include <vector>
#include <array>

namespace crispr_gpu {

float score_mismatch_count(uint8_t mismatches);
float score_mit(const std::vector<uint8_t> &mismatch_positions, uint8_t guide_length);
float score_cfd_bits(uint64_t guide_bits, uint64_t site_bits, uint8_t guide_length);
void load_cfd_tables(const std::string &json_path);
void ensure_default_tables_loaded();

struct ScoringTables {
  std::array<float,20> mit_position_penalty{};
  std::array<float,20> cfd_position_weight{};
  float cfd_type_weight[4][4]{};
};

// Returns process-wide tables, loading defaults or a custom path once per process.
const ScoringTables &get_scoring_tables(const ScoreParams &params);

} // namespace crispr_gpu
