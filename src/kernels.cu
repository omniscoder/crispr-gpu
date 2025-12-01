#include <cuda_runtime.h>
#include <algorithm>
#include <cstring>
#include "crispr_gpu/types.hpp"
#include "crispr_gpu/scoring.hpp"
#include "crispr_gpu/engine.hpp"

namespace crispr_gpu {
namespace detail {

struct DeviceScoringTables {
  float mit_position_penalty[20];
  float cfd_position_weight[20];
  float cfd_type_weight[4][4];
};

__device__ __constant__ DeviceScoringTables kDeviceTables;

__device__ inline uint8_t mismatch_count_2bit(uint64_t guide_bits, uint64_t site_bits, uint8_t guide_length) {
  uint8_t used_bits = static_cast<uint8_t>(guide_length * 2);
  uint64_t mask = (used_bits == 64) ? ~0ULL : ((1ULL << used_bits) - 1ULL);
  uint64_t x = (guide_bits ^ site_bits) & mask;
  uint64_t hi = x & 0xAAAAAAAAAAAAAAAAull;
  uint64_t lo = x & 0x5555555555555555ull;
  uint64_t collapsed = (hi >> 1) | lo;
  uint64_t mm_mask = 0x5555555555555555ull;
  uint64_t mm = collapsed & mm_mask;
  return static_cast<uint8_t>(__popcll(mm));
}

__device__ inline float score_hamming(uint8_t mismatches) {
  return 1.0f / (1.0f + static_cast<float>(mismatches));
}

__device__ inline float score_mit_device(const uint8_t *pos, uint8_t count) {
  float score = 1.0f;
  for (uint8_t i = 0; i < count; ++i) {
    uint8_t p = pos[i];
    float pen = (p < 20) ? kDeviceTables.mit_position_penalty[p] : 0.5f;
    score *= (1.0f - pen);
    if (i > 0 && pos[i] - pos[i - 1] <= 1) score *= 0.8f;
  }
  return score;
}

__device__ inline float score_cfd_device(uint64_t guide_bits,
                                         uint64_t site_bits,
                                         uint8_t guide_length) {
  float score = 1.0f;
  for (uint8_t p = 0; p < guide_length; ++p) {
    uint8_t shift = static_cast<uint8_t>(2 * (guide_length - 1 - p));
    uint8_t g = static_cast<uint8_t>((guide_bits >> shift) & 0b11);
    uint8_t s = static_cast<uint8_t>((site_bits >> shift) & 0b11);
    if (g == s) continue;
    float pos_w = (p < 20) ? kDeviceTables.cfd_position_weight[p] : 0.7f;
    float type_w = kDeviceTables.cfd_type_weight[g][s];
    score *= (pos_w * type_w);
  }
  return score;
}

__global__ void off_target_kernel(const SiteRecord *__restrict__ sites,
                                  uint32_t num_sites,
                                  uint64_t guide_bits,
                                  uint8_t max_mm,
                                  uint8_t guide_length,
                                  ScoreModel score_model,
                                  DeviceHit *__restrict__ out_hits,
                                  uint32_t *__restrict__ out_count) {
  uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for (uint32_t i = tid; i < num_sites; i += stride) {
    const SiteRecord site = sites[i];
    uint8_t mm = mismatch_count_2bit(guide_bits, site.seq_bits, guide_length);
    if (mm > max_mm) continue;

    float score = 0.0f;
    if (score_model == ScoreModel::Hamming) {
      score = score_hamming(mm);
    } else if (score_model == ScoreModel::MIT) {
      uint8_t posbuf[32];
      uint8_t mcount = 0;
      uint64_t a = guide_bits;
      uint64_t b = site.seq_bits;
      for (uint8_t p = 0; p < guide_length; ++p) {
        uint8_t shift = static_cast<uint8_t>(2 * (guide_length - 1 - p));
        uint8_t aa = static_cast<uint8_t>((a >> shift) & 0b11);
        uint8_t bb = static_cast<uint8_t>((b >> shift) & 0b11);
        if (aa != bb) posbuf[mcount++] = p;
      }
      score = score_mit_device(posbuf, mcount);
    } else { // CFD
      score = score_cfd_device(guide_bits, site.seq_bits, guide_length);
    }
    uint32_t idx = atomicAdd(out_count, 1u);
    out_hits[idx].site_index = i;
    out_hits[idx].mismatches = mm;
    out_hits[idx].score = score;
  }
}

void upload_scoring_tables(const ScoringTables &tables) {
  static DeviceScoringTables cached{};
  static bool initialised = false;
  DeviceScoringTables tmp{};
  std::memcpy(tmp.mit_position_penalty, tables.mit_position_penalty.data(), sizeof(tmp.mit_position_penalty));
  std::memcpy(tmp.cfd_position_weight, tables.cfd_position_weight.data(), sizeof(tmp.cfd_position_weight));
  std::memcpy(tmp.cfd_type_weight, tables.cfd_type_weight, sizeof(tmp.cfd_type_weight));

  if (!initialised || std::memcmp(&cached, &tmp, sizeof(DeviceScoringTables)) != 0) {
    cached = tmp;
    cudaMemcpyToSymbol(kDeviceTables, &cached, sizeof(DeviceScoringTables));
    initialised = true;
  }
}

void launch_off_target_kernel(const SiteRecord *d_sites,
                              uint32_t num_sites,
                              uint64_t guide_bits,
                              uint8_t max_mm,
                              uint8_t guide_length,
                              ScoreModel score_model,
                              DeviceHit *d_hits,
                              uint32_t *d_count,
                              cudaStream_t stream) {
  int block = 256;
  int grid = std::min<int>((num_sites + block - 1) / block, 65535);
  off_target_kernel<<<grid, block, 0, stream>>>(d_sites, num_sites, guide_bits, max_mm,
                                                guide_length, score_model, d_hits, d_count);
}

} // namespace detail
} // namespace crispr_gpu
