#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <chrono>
#include <algorithm>

// Simple device-side mismatch count for 2-bit encoded bases (20 bp assumed).
__device__ __forceinline__ uint8_t mismatch_count_2bit(uint64_t guide_bits, uint64_t site_bits) {
  uint64_t x = guide_bits ^ site_bits;
  uint64_t mism = (x | (x >> 1)) & 0x5555555555555555ULL;
  return static_cast<uint8_t>(__popcll(mism));
}

__global__ void hamming_kernel(const uint64_t *sites, uint64_t guide_bits,
                               uint8_t guide_length, uint64_t iters,
                               uint64_t N, unsigned long long *out_sum) {
  uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  uint64_t stride = static_cast<uint64_t>(blockDim.x) * gridDim.x;
  uint64_t acc = 0;
  for (uint64_t iter = 0; iter < iters; ++iter) {
    for (uint64_t i = tid; i < N; i += stride) {
      uint8_t mm = mismatch_count_2bit(guide_bits, sites[i]);
      acc += mm;
    }
  }
  atomicAdd(out_sum, static_cast<unsigned long long>(acc));
}

int main() {
  const uint64_t default_N = 10'000'000ULL;
  const uint64_t default_iters = 50;
  uint64_t N = std::strtoull(std::getenv("MB_N") ? std::getenv("MB_N") : "0", nullptr, 10);
  uint64_t iters = std::strtoull(std::getenv("MB_ITERS") ? std::getenv("MB_ITERS") : "0", nullptr, 10);
  if (N == 0) N = default_N;
  if (iters == 0) iters = default_iters;

  uint64_t guide_bits = 0xAAAAAAAAAAAAAAAULL; // arbitrary pattern

  uint64_t *d_sites = nullptr;
  unsigned long long *d_sum = nullptr;

  cudaMalloc(&d_sites, N * sizeof(uint64_t));
  cudaMalloc(&d_sum, sizeof(unsigned long long));

  // Fill sites with pseudo data (here zeros is fine)
  cudaMemset(d_sites, 0, N * sizeof(uint64_t));
  unsigned long long zero = 0;
  cudaMemcpy(d_sum, &zero, sizeof(unsigned long long), cudaMemcpyHostToDevice);

  dim3 block(256);
  uint64_t grid_x = (N + block.x - 1) / block.x;
  if (grid_x > 65535ULL) grid_x = 65535ULL;
  dim3 grid(static_cast<unsigned int>(grid_x));

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);
  hamming_kernel<<<grid, block>>>(d_sites, guide_bits, 20, iters, N, d_sum);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);

  float ms = 0.0f;
  cudaEventElapsedTime(&ms, start, stop);

  unsigned long long acc_host = 0;
  cudaMemcpy(&acc_host, d_sum, sizeof(unsigned long long), cudaMemcpyDeviceToHost);

  double seconds = ms / 1000.0;
  double candidates = static_cast<double>(N) * static_cast<double>(iters);
  double cps = candidates / seconds;

  std::printf("kernel_candidates_per_sec: %.3f\n", cps);
  std::printf("total_candidates: %.0f\n", candidates);
  std::printf("time_sec: %.6f\n", seconds);
  std::printf("accum: %llu\n", static_cast<unsigned long long>(acc_host));

  cudaFree(d_sites);
  cudaFree(d_sum);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  return 0;
}
