#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <stdexcept>
#include <array>

namespace crispr_gpu {

enum class ScoreModel {
  Hamming,
  MIT,
  CFD
};

enum class Backend {
  CPU,
  GPU
};

// Stage-1 search backend: today brute-force scan; FM-index/seeded search planned.
enum class SearchBackend {
  BruteForce,
  FMIndex
};

struct ScoreParams {
  ScoreModel model{ScoreModel::Hamming};
  std::string table_path{}; // optional override for MIT/CFD tables
};

struct EngineParams {
  ScoreParams score_params{};
  uint8_t max_mismatches{4};
  Backend backend{
#ifdef CRISPR_GPU_ENABLE_CUDA
      Backend::GPU
#else
      Backend::CPU
#endif
  };
  SearchBackend search_backend{SearchBackend::BruteForce};
};

struct Guide {
  std::string name;
  std::string sequence;  // raw bases
  std::string pam{"NGG"};
};

struct EncodedGuide {
  uint64_t bits{0};
  uint8_t length{0};
  std::string name;
};

struct alignas(16) SiteRecord {
  uint64_t seq_bits{0};
  uint32_t chrom_id{0};
  uint32_t pos{0};
  uint8_t strand{0}; // 0 = '+', 1 = '-'
  uint8_t pad[3]{0, 0, 0};
};

struct ChromInfo {
  std::string name;
  uint64_t length{0};
};

struct IndexMeta {
  uint8_t guide_length{20};
  std::string pam{"NGG"};
  bool both_strands{true};
};

struct OffTargetHit {
  std::string guide_name;
  uint32_t chrom_id{0};
  uint32_t pos{0};
  char strand{'+'};
  uint8_t mismatches{0};
  float score{0.0f};
};

struct BenchStat {
  std::uint64_t candidates{0};
  double seconds{0.0};
  double cgct() const { return seconds > 0.0 ? static_cast<double>(candidates) / seconds : 0.0; }
};

// Utility helpers
uint8_t base_to_bits(char b);
char bits_to_base(uint8_t bits);
uint64_t encode_sequence_2bit(const std::string &seq);
std::string decode_sequence_2bit(uint64_t bits, uint8_t length);
std::string revcomp(const std::string &seq);

uint8_t hamming_distance_2bit(uint64_t a, uint64_t b, uint8_t length);

} // namespace crispr_gpu
