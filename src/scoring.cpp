#include "crispr_gpu/scoring.hpp"

#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include <cstdlib>
#include <cstdint>

namespace crispr_gpu {

float score_mismatch_count(uint8_t mismatches) {
  return 1.0f / (1.0f + static_cast<float>(mismatches));
}

static ScoringTables g_tables = {
    /*mit_position_penalty*/ {0.1f, 0.1f, 0.1f, 0.12f, 0.12f, 0.14f, 0.14f, 0.16f, 0.16f, 0.18f,
                              0.18f, 0.20f, 0.22f, 0.24f, 0.26f, 0.30f, 0.34f, 0.38f, 0.42f, 0.45f},
    /*cfd_position_weight*/  {0.61f, 0.61f, 0.62f, 0.63f, 0.64f, 0.65f, 0.66f, 0.67f, 0.68f, 0.69f,
                              0.70f, 0.71f, 0.72f, 0.74f, 0.76f, 0.78f, 0.80f, 0.83f, 0.86f, 0.90f},
    /*cfd_type_weight*/      {{1.0f, 0.97f, 0.48f, 0.88f},
                              {0.87f, 1.0f, 0.55f, 0.73f},
                              {0.14f, 0.34f, 1.0f, 0.03f},
                              {0.76f, 0.67f, 0.32f, 1.0f}}
};

void load_cfd_tables(const std::string &json_path) {
  std::ifstream in(json_path);
  if (!in) throw std::runtime_error("Failed to open CFD table: " + json_path);
  nlohmann::json j;
  in >> j;
  if (j.contains("mm_scores")) {
    auto mm = j["mm_scores"];
    for (auto &kv : mm.items()) {
      // key format "AC": guide A, genome C
      const std::string &k = kv.key();
      if (k.size() != 2) continue;
      int g = std::string("ACGT").find(k[0]);
      int s = std::string("ACGT").find(k[1]);
      if (g >= 0 && s >= 0) {
        g_tables.cfd_type_weight[g][s] = kv.value().get<float>();
      }
    }
  }
  if (j.contains("position")) {
    auto pos = j["position"];
    for (size_t i = 0; i < std::min<size_t>(g_tables.cfd_position_weight.size(), pos.size()); ++i) {
      g_tables.cfd_position_weight[i] = pos[i].get<float>();
    }
  }
  if (j.contains("mit_position_penalty")) {
    auto pos = j["mit_position_penalty"];
    for (size_t i = 0; i < std::min<size_t>(g_tables.mit_position_penalty.size(), pos.size()); ++i) {
      g_tables.mit_position_penalty[i] = pos[i].get<float>();
    }
  }
}

static std::string default_cfd_path() {
#ifdef CRISPR_GPU_DATA_DIR
  return std::string(CRISPR_GPU_DATA_DIR) + "/cfd_default.json";
#else
  return "";
#endif
}

void ensure_default_tables_loaded() {
  static bool loaded = false;
  if (loaded) return;
  std::string path = default_cfd_path();
  if (!path.empty()) {
    try {
      load_cfd_tables(path);
    } catch (...) {
      // keep defaults
    }
  }
  loaded = true;
}

float score_mit(const std::vector<uint8_t> &mismatch_positions, [[maybe_unused]] uint8_t guide_length) {
  float score = 1.0f;
  for (auto pos : mismatch_positions) {
    if (pos < g_tables.mit_position_penalty.size()) {
      score *= (1.0f - g_tables.mit_position_penalty[pos]);
    } else {
      score *= 0.5f; // conservative penalty for unexpected length
    }
  }
  // Additional penalty for mismatch clustering (crude approximation)
  for (size_t i = 1; i < mismatch_positions.size(); ++i) {
    if (mismatch_positions[i] - mismatch_positions[i - 1] <= 1) {
      score *= 0.8f;
    }
  }
  return score;
}

static inline uint8_t bits_at(uint64_t bits, uint8_t pos_from_msb, uint8_t guide_length) {
  uint8_t shift = static_cast<uint8_t>(2 * (guide_length - 1 - pos_from_msb));
  return static_cast<uint8_t>((bits >> shift) & 0b11);
}

float score_cfd_bits(uint64_t guide_bits, uint64_t site_bits, uint8_t guide_length) {
  float score = 1.0f;
  for (uint8_t p = 0; p < guide_length; ++p) {
    uint8_t g = bits_at(guide_bits, p, guide_length);
    uint8_t s = bits_at(site_bits, p, guide_length);
    if (g == s) continue;
    float pos_w = (p < g_tables.cfd_position_weight.size()) ? g_tables.cfd_position_weight[p] : 0.7f;
    float type_w = g_tables.cfd_type_weight[g][s];
    score *= (pos_w * type_w);
  }
  return score;
}

const ScoringTables &get_scoring_tables(const ScoreParams &params) {
  static std::string last_path;
  if (!params.table_path.empty() && params.table_path != last_path) {
    load_cfd_tables(params.table_path);
    last_path = params.table_path;
    return g_tables;
  }
  ensure_default_tables_loaded();
  return g_tables;
}

} // namespace crispr_gpu
