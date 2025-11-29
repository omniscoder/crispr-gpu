#pragma once

#include "crispr_gpu/types.hpp"
#include <vector>

namespace crispr_gpu {

float score_mismatch_count(uint8_t mismatches);
float score_mit(const std::vector<uint8_t> &mismatch_positions, uint8_t guide_length);
float score_cfd_bits(uint64_t guide_bits, uint64_t site_bits, uint8_t guide_length);
void load_cfd_tables(const std::string &json_path);
void ensure_default_tables_loaded();

} // namespace crispr_gpu
