#include <cuda_runtime.h>
#include <algorithm>
#include "crispr_gpu/types.hpp"
#include "crispr_gpu/scoring.hpp"
#include "crispr_gpu/engine.hpp"

namespace crispr_gpu {
namespace detail {

__device__ __constant__ float kMitPositionPenalty[20] = {
    0.1f, 0.1f, 0.1f, 0.12f, 0.12f, 0.14f, 0.14f, 0.16f, 0.16f, 0.18f,
    0.18f, 0.20f, 0.22f, 0.24f, 0.26f, 0.30f, 0.34f, 0.38f, 0.42f, 0.45f};

__device__ inline uint8_t mismatch_count_2bit(uint64_t guide_bits, uint64_t site_bits) {
  uint64_t x = guide_bits ^ site_bits;
  uint64_t hi = x & 0xAAAAAAAAAAAAAAAAull;
  uint64_t lo = x & 0x5555555555555555ull;
  uint64_t collapsed = (hi >> 1) | lo;
  uint64_t mask = 0x5555555555555555ull;
  uint64_t mm = collapsed & mask;
  return static_cast<uint8_t>(__popcll(mm));
}

__device__ inline float score_hamming(uint8_t mismatches) {
  return 1.0f / (1.0f + static_cast<float>(mismatches));
}

__device__ inline float score_mit_device(const uint8_t *pos, uint8_t count) {
  float score = 1.0f;
  for (uint8_t i = 0; i < count; ++i) {
    uint8_t p = pos[i];
    float pen = (p < 20) ? kMitPositionPenalty[p] : 0.5f;
    score *= (1.0f - pen);
    if (i > 0 && pos[i] - pos[i - 1] <= 1) score *= 0.8f;
  }
  return score;
}

__device__ inline float score_cfd_device(const uint8_t *pos, uint8_t count) {
  float score = 1.0f;
  for (uint8_t i = 0; i < count; ++i) {
    float weight = 0.8f - 0.02f * static_cast<float>(pos[i]);
    if (weight < 0.4f) weight = 0.4f;
    score *= weight;
  }
  return score;
}

__device__ inline float score_hit_device(const uint8_t *mm_positions, uint8_t mm_count,
                                         ScoreParams score_params) {
  switch (score_params.model) {
    case ScoreModel::Hamming:
      return score_hamming(mm_count);
    case ScoreModel::MIT:
      return score_mit_device(mm_positions, mm_count);
    case ScoreModel::CFD:
      return score_cfd_device(mm_positions, mm_count);
  }
  return 0.0f;
}

__global__ void off_target_kernel(const SiteRecord *sites,
                                  uint32_t num_sites,
                                  uint64_t guide_bits,
                                  uint8_t max_mm,
                                  uint8_t guide_length,
                                  ScoreParams score_params,
                                  DeviceHit *out_hits,
                                  uint32_t *out_count) {
  uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for (uint32_t i = tid; i < num_sites; i += stride) {
    const SiteRecord site = sites[i];
    uint8_t mm = mismatch_count_2bit(guide_bits, site.seq_bits);
    if (mm > max_mm) continue;

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

    float score = score_hit_device(posbuf, mcount, score_params);
    uint32_t idx = atomicAdd(out_count, 1u);
    out_hits[idx].site_index = i;
    out_hits[idx].mismatches = mm;
    out_hits[idx].score = score;
  }
}

void launch_off_target_kernel(const SiteRecord *d_sites,
                              uint32_t num_sites,
                              uint64_t guide_bits,
                              uint8_t max_mm,
                              uint8_t guide_length,
                              ScoreParams score_params,
                              DeviceHit *d_hits,
                              uint32_t *d_count,
                              cudaStream_t stream) {
  int block = 256;
  int grid = std::min<int>((num_sites + block - 1) / block, 65535);
  off_target_kernel<<<grid, block, 0, stream>>>(d_sites, num_sites, guide_bits, max_mm,
                                                guide_length, score_params, d_hits, d_count);
}

} // namespace detail
} // namespace crispr_gpu
