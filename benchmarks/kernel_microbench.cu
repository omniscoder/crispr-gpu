#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>
#include <chrono>

__device__ __forceinline__ uint8_t mismatch_count_2bit(uint64_t guide_bits, uint64_t site_bits) {
  uint64_t x = guide_bits ^ site_bits;
  uint64_t mism = (x | (x >> 1)) & 0x5555555555555555ULL;
  return static_cast<uint8_t>(__popcll(mism));
}

struct DeviceHit { uint32_t site_index; uint8_t mismatches; float score; };

__global__ void off_target_kernel(const uint64_t *sites,
                                  uint32_t num_sites,
                                  uint64_t guide_bits,
                                  uint8_t max_mm,
                                  DeviceHit *out_hits,
                                  uint32_t *out_count) {
  uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for (uint32_t i = tid; i < num_sites; i += stride) {
    uint64_t site_bits = sites[i];
    uint8_t mm = mismatch_count_2bit(guide_bits, site_bits);
    if (mm > max_mm) continue;
    uint32_t idx = atomicAdd(out_count, 1u);
    out_hits[idx].site_index = i;
    out_hits[idx].mismatches = mm;
    out_hits[idx].score = 1.0f; // dummy
  }
}

int main() {
  const uint32_t num_sites = 6'246'000; // matches bench index size
  const uint64_t guide_bits = 0x1b1b1b1b1bULL; // arbitrary 20bp pattern
  const int block = 256;
  const int grid = (num_sites + block - 1) / block; // full occupancy pattern used in prod

  uint64_t *d_sites = nullptr;
  DeviceHit *d_hits = nullptr;
  uint32_t *d_count = nullptr;

  cudaMalloc(&d_sites, num_sites * sizeof(uint64_t));
  cudaMalloc(&d_hits, num_sites * sizeof(DeviceHit)); // upper bound
  cudaMalloc(&d_count, sizeof(uint32_t));
  cudaMemset(d_sites, 0, num_sites * sizeof(uint64_t));
  cudaMemset(d_count, 0, sizeof(uint32_t));

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);
  off_target_kernel<<<grid, block>>>(d_sites, num_sites, guide_bits, /*max_mm=*/4, d_hits, d_count);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);

  float ms = 0.0f;
  cudaEventElapsedTime(&ms, start, stop);

  uint32_t count_host = 0;
  cudaMemcpy(&count_host, d_count, sizeof(uint32_t), cudaMemcpyDeviceToHost);

  double seconds = ms / 1000.0;
  double candidates = static_cast<double>(num_sites);
  double cgct = candidates / seconds;

  std::printf("grid %d block %d\n", grid, block);
  std::printf("time_sec %.6f\n", seconds);
  std::printf("total_candidates %.0f\n", candidates);
  std::printf("cgct_candidates_per_sec %.3e\n", cgct);
  std::printf("hits %u\n", count_host);

  cudaFree(d_sites);
  cudaFree(d_hits);
  cudaFree(d_count);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  return 0;
}
