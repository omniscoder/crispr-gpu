#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

#define CUDA_CHECK(call)                                                                              \
  do {                                                                                                \
    cudaError_t _err = (call);                                                                         \
    if (_err != cudaSuccess) {                                                                         \
      std::fprintf(stderr, "CUDA error %s (%d) at %s:%d: %s\n", cudaGetErrorName(_err),                \
                   static_cast<int>(_err), __FILE__, __LINE__, cudaGetErrorString(_err));             \
      std::exit(1);                                                                                   \
    }                                                                                                 \
  } while (0)

__device__ __forceinline__ uint8_t mismatch_count_2bit(uint64_t guide_bits, uint64_t site_bits) {
  uint64_t x = guide_bits ^ site_bits;
  uint64_t mism = (x | (x >> 1)) & 0x5555555555555555ULL;
  return static_cast<uint8_t>(__popcll(mism));
}

struct DeviceHit { uint32_t site_index; uint8_t mismatches; float score; };

__global__ void init_sites_kernel(uint64_t *sites,
                                  uint32_t num_sites,
                                  uint64_t guide_bits,
                                  uint64_t nonhit_bits,
                                  uint32_t hit_sites) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_sites) return;
  sites[i] = (i < hit_sites) ? guide_bits : nonhit_bits;
}

__global__ void off_target_kernel(const uint64_t *sites,
                                  uint32_t num_sites,
                                  uint64_t guide_bits,
                                  uint8_t max_mm,
                                  DeviceHit *out_hits,
                                  uint32_t hit_capacity,
                                  uint32_t *out_count,
                                  uint32_t *out_overflow) {
  uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  uint32_t stride = blockDim.x * gridDim.x;

  for (uint32_t i = tid; i < num_sites; i += stride) {
    uint64_t site_bits = sites[i];
    uint8_t mm = mismatch_count_2bit(guide_bits, site_bits);
    if (mm > max_mm) continue;
    uint32_t idx = atomicAdd(out_count, 1u);
    if (idx < hit_capacity) {
      out_hits[idx].site_index = i;
      out_hits[idx].mismatches = mm;
      out_hits[idx].score = 1.0f; // dummy
    } else {
      atomicAdd(out_overflow, 1u);
    }
  }
}

static uint32_t parse_u32(const char *s, const char *name) {
  char *end = nullptr;
  unsigned long v = std::strtoul(s, &end, 10);
  if (!s[0] || (end && *end) || v > 0xFFFFFFFFu) {
    std::fprintf(stderr, "Invalid %s: %s\n", name, s);
    std::exit(2);
  }
  return static_cast<uint32_t>(v);
}

static double parse_f64(const char *s, const char *name) {
  char *end = nullptr;
  double v = std::strtod(s, &end);
  if (!s[0] || (end && *end)) {
    std::fprintf(stderr, "Invalid %s: %s\n", name, s);
    std::exit(2);
  }
  return v;
}

static void print_usage() {
  std::fprintf(stderr,
               "kernel_microbench\n"
               "  Device-only kernel throughput benchmark (single-guide scan).\n\n"
               "Usage:\n"
               "  kernel_microbench [--num-sites N] [--block N] [--grid N] [--max-mm K]\n"
               "                  [--hit-fraction X] [--iters N] [--warmup N]\n"
               "                  [--hit-capacity N] [--device ID]\n"
               "                  [--format kv|json] [--output PATH]\n\n"
               "Notes:\n"
               "  --hit-fraction controls how many sites are exact matches (mm=0).\n"
               "  Non-hit sites are constructed to have mm=20 for a 20nt guide.\n");
}

int main(int argc, char **argv) {
  uint32_t num_sites = 6'246'000;              // matches 50 Mb synthetic index size
  uint8_t max_mm = 4;
  uint32_t block = 256;
  uint32_t grid = 0;                           // auto by default
  uint32_t iters = 20;
  uint32_t warmup = 2;
  double hit_fraction = 0.0;                   // 0 => no atomics/writes, pure scan + popc
  uint32_t hit_capacity = 0;                   // 0 => num_sites
  int device = 0;
  std::string format = "kv";
  std::string output_path;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      print_usage();
      return 0;
    }
    auto need_value = [&](const char *name) -> const char * {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "Missing value for %s\n", name);
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "--num-sites") { num_sites = parse_u32(need_value("--num-sites"), "--num-sites"); continue; }
    if (arg == "--block") { block = parse_u32(need_value("--block"), "--block"); continue; }
    if (arg == "--grid") { grid = parse_u32(need_value("--grid"), "--grid"); continue; }
    if (arg == "--max-mm") { max_mm = static_cast<uint8_t>(parse_u32(need_value("--max-mm"), "--max-mm")); continue; }
    if (arg == "--iters") { iters = parse_u32(need_value("--iters"), "--iters"); continue; }
    if (arg == "--warmup") { warmup = parse_u32(need_value("--warmup"), "--warmup"); continue; }
    if (arg == "--hit-capacity") { hit_capacity = parse_u32(need_value("--hit-capacity"), "--hit-capacity"); continue; }
    if (arg == "--hit-fraction") { hit_fraction = parse_f64(need_value("--hit-fraction"), "--hit-fraction"); continue; }
    if (arg == "--device") { device = static_cast<int>(parse_u32(need_value("--device"), "--device")); continue; }
    if (arg == "--format") { format = need_value("--format"); continue; }
    if (arg == "--output") { output_path = need_value("--output"); continue; }
    std::fprintf(stderr, "Unknown argument: %s\n", arg.c_str());
    print_usage();
    return 2;
  }

  if (hit_fraction < 0.0 || hit_fraction > 1.0) {
    std::fprintf(stderr, "--hit-fraction must be within [0,1]\n");
    return 2;
  }
  if (block == 0) {
    std::fprintf(stderr, "--block must be > 0\n");
    return 2;
  }
  if (iters == 0) {
    std::fprintf(stderr, "--iters must be > 0\n");
    return 2;
  }

  CUDA_CHECK(cudaSetDevice(device));

  cudaDeviceProp prop{};
  CUDA_CHECK(cudaGetDeviceProperties(&prop, device));

  int runtime_version = 0;
  int driver_version = 0;
  CUDA_CHECK(cudaRuntimeGetVersion(&runtime_version));
  CUDA_CHECK(cudaDriverGetVersion(&driver_version));

  const uint8_t guide_len = 20;
  const uint64_t guide_bits = 0x1b1b1b1b1bULL; // arbitrary 20bp pattern (packed 2-bit)
  const uint64_t guide_mask = (1ULL << (2 * guide_len)) - 1ULL;
  const uint64_t nonhit_bits = (guide_bits ^ guide_mask);
  const uint32_t hit_sites = static_cast<uint32_t>(hit_fraction * static_cast<double>(num_sites));
  if (hit_capacity == 0) hit_capacity = num_sites;

  if (grid == 0) {
    grid = (num_sites + block - 1) / block;
    if (grid == 0) grid = 1;
  }

  uint64_t *d_sites = nullptr;
  DeviceHit *d_hits = nullptr;
  uint32_t *d_count = nullptr;
  uint32_t *d_overflow = nullptr;

  CUDA_CHECK(cudaMalloc(&d_sites, static_cast<size_t>(num_sites) * sizeof(uint64_t)));
  CUDA_CHECK(cudaMalloc(&d_hits, static_cast<size_t>(hit_capacity) * sizeof(DeviceHit)));
  CUDA_CHECK(cudaMalloc(&d_count, sizeof(uint32_t)));
  CUDA_CHECK(cudaMalloc(&d_overflow, sizeof(uint32_t)));

  {
    const uint32_t init_block = 256;
    const uint32_t init_grid = (num_sites + init_block - 1) / init_block;
    init_sites_kernel<<<init_grid, init_block>>>(d_sites, num_sites, guide_bits, nonhit_bits, hit_sites);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  cudaEvent_t start{}, stop{};
  CUDA_CHECK(cudaEventCreate(&start));
  CUDA_CHECK(cudaEventCreate(&stop));

  std::vector<float> times_ms;
  times_ms.reserve(iters);

  for (uint32_t i = 0; i < warmup + iters; ++i) {
    CUDA_CHECK(cudaMemset(d_count, 0, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemset(d_overflow, 0, sizeof(uint32_t)));

    CUDA_CHECK(cudaEventRecord(start));
    off_target_kernel<<<grid, block>>>(d_sites, num_sites, guide_bits, max_mm, d_hits, hit_capacity, d_count, d_overflow);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));

    float ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&ms, start, stop));
    if (i >= warmup) times_ms.push_back(ms);
  }

  uint32_t count_host = 0;
  uint32_t overflow_host = 0;
  CUDA_CHECK(cudaMemcpy(&count_host, d_count, sizeof(uint32_t), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(&overflow_host, d_overflow, sizeof(uint32_t), cudaMemcpyDeviceToHost));

  std::vector<float> sorted = times_ms;
  std::sort(sorted.begin(), sorted.end());
  const float min_ms = sorted.front();
  const float max_ms = sorted.back();
  const float median_ms = sorted[sorted.size() / 2];
  const double mean_ms = std::accumulate(times_ms.begin(), times_ms.end(), 0.0) / static_cast<double>(times_ms.size());

  const double mean_sec = mean_ms / 1000.0;
  const double candidates = static_cast<double>(num_sites);
  const double cgct = candidates / mean_sec;

  FILE *out = stdout;
  FILE *owned = nullptr;
  if (!output_path.empty()) {
    owned = std::fopen(output_path.c_str(), "wb");
    if (!owned) {
      std::perror("fopen");
      return 1;
    }
    out = owned;
  }

  if (format == "json") {
    std::fprintf(out,
                 "{\n"
                 "  \"schema\": \"crispr-gpu/kernel_microbench/v1\",\n"
                 "  \"schema_version\": 1,\n"
                 "  \"device\": {\n"
                 "    \"id\": %d,\n"
                 "    \"name\": \"%s\",\n"
                 "    \"cc\": \"%d.%d\",\n"
                 "    \"sm_count\": %d,\n"
                 "    \"clock_khz\": %d,\n"
                 "    \"driver_version\": %d,\n"
                 "    \"runtime_version\": %d\n"
                 "  },\n"
                 "  \"kernel\": {\n"
                 "    \"num_sites\": %u,\n"
                 "    \"guide_length\": %u,\n"
                 "    \"max_mm\": %u,\n"
                 "    \"block\": %u,\n"
                 "    \"grid\": %u,\n"
                 "    \"hit_fraction\": %.6f,\n"
                 "    \"hit_capacity\": %u\n"
                 "  },\n"
                 "  \"timing\": {\n"
                 "    \"iters\": %u,\n"
                 "    \"warmup\": %u,\n"
                 "    \"min_ms\": %.6f,\n"
                 "    \"median_ms\": %.6f,\n"
                 "    \"mean_ms\": %.6f,\n"
                 "    \"max_ms\": %.6f\n"
                 "  },\n"
                 "  \"results\": {\n"
                 "    \"mean_sec\": %.9f,\n"
                 "    \"cgct_candidates_per_sec\": %.6e,\n"
                 "    \"hits\": %u,\n"
                 "    \"overflow\": %u\n"
                 "  }\n"
                 "}\n",
                 device, prop.name, prop.major, prop.minor, prop.multiProcessorCount,
                 prop.clockRate, driver_version, runtime_version, num_sites,
                 static_cast<unsigned>(guide_len), static_cast<unsigned>(max_mm), block, grid,
                 hit_fraction, hit_capacity, iters, warmup,
                 static_cast<double>(min_ms), static_cast<double>(median_ms), mean_ms, static_cast<double>(max_ms),
                 mean_sec, cgct, count_host, overflow_host);
  } else if (format == "kv") {
    std::fprintf(out, "device_id %d\n", device);
    std::fprintf(out, "device_name %s\n", prop.name);
    std::fprintf(out, "device_cc %d.%d\n", prop.major, prop.minor);
    std::fprintf(out, "device_sm_count %d\n", prop.multiProcessorCount);
    std::fprintf(out, "cuda_driver_version %d\n", driver_version);
    std::fprintf(out, "cuda_runtime_version %d\n", runtime_version);
    std::fprintf(out, "num_sites %u\n", num_sites);
    std::fprintf(out, "guide_length %u\n", static_cast<unsigned>(guide_len));
    std::fprintf(out, "max_mm %u\n", static_cast<unsigned>(max_mm));
    std::fprintf(out, "grid %u\n", grid);
    std::fprintf(out, "block %u\n", block);
    std::fprintf(out, "hit_fraction %.6f\n", hit_fraction);
    std::fprintf(out, "hit_capacity %u\n", hit_capacity);
    std::fprintf(out, "iters %u\n", iters);
    std::fprintf(out, "warmup %u\n", warmup);
    std::fprintf(out, "min_ms %.6f\n", static_cast<double>(min_ms));
    std::fprintf(out, "median_ms %.6f\n", static_cast<double>(median_ms));
    std::fprintf(out, "mean_ms %.6f\n", mean_ms);
    std::fprintf(out, "max_ms %.6f\n", static_cast<double>(max_ms));
    std::fprintf(out, "mean_sec %.9f\n", mean_sec);
    std::fprintf(out, "cgct_candidates_per_sec %.6e\n", cgct);
    std::fprintf(out, "hits %u\n", count_host);
    std::fprintf(out, "overflow %u\n", overflow_host);
  } else {
    std::fprintf(stderr, "Unknown --format: %s (expected kv|json)\n", format.c_str());
    if (owned) std::fclose(owned);
    return 2;
  }

  if (owned) std::fclose(owned);

  CUDA_CHECK(cudaFree(d_sites));
  CUDA_CHECK(cudaFree(d_hits));
  CUDA_CHECK(cudaFree(d_count));
  CUDA_CHECK(cudaFree(d_overflow));
  CUDA_CHECK(cudaEventDestroy(start));
  CUDA_CHECK(cudaEventDestroy(stop));
  return 0;
}
