#include "crispr_gpu/engine.hpp"
#include "crispr_gpu/types.hpp"
#include "crispr_gpu/scoring.hpp"

#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cstdio>

#ifdef CRISPR_GPU_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

namespace crispr_gpu {

#ifdef CRISPR_GPU_ENABLE_CUDA
namespace detail {
void launch_off_target_kernel(const SiteRecord *d_sites,
                              uint32_t num_sites,
                              uint64_t guide_bits,
                              uint8_t max_mm,
                              uint8_t guide_length,
                              ScoreParams score_params,
                              DeviceHit *d_hits,
                              uint32_t *d_count,
                              cudaStream_t stream);
}                                                                                
#endif

namespace {

bool timing_enabled() {
  static int cached = -1;
  if (cached == -1) {
    const char *env = std::getenv("CRISPR_GPU_TIMING");
    cached = (env && env[0] != '0') ? 1 : 0;
  }
  return cached == 1;
}

struct ScopedTimer {
  const char *name;
  bool active;
  std::chrono::steady_clock::time_point start;
  ScopedTimer(const char *n, bool on) : name(n), active(on), start(std::chrono::steady_clock::now()) {}
  ~ScopedTimer() {
    if (!active) return;
    auto end = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::fprintf(stderr, "[timing] %s: %.3f ms\n", name, ms);
  }
};

EncodedGuide encode_guide(const Guide &guide, uint8_t expected_length) {
  if (guide.sequence.size() != expected_length) {
    throw std::runtime_error("Guide length mismatch vs index");
  }
  EncodedGuide eg;
  eg.bits = encode_sequence_2bit(guide.sequence);
  eg.length = expected_length;
  eg.name = guide.name;
  return eg;
}

std::vector<uint8_t> mismatch_positions(uint64_t a, uint64_t b, uint8_t length) {
  std::vector<uint8_t> positions;
  for (uint8_t i = 0; i < length; ++i) {
    uint8_t shift = static_cast<uint8_t>(2 * (length - 1 - i));
    uint8_t aa = static_cast<uint8_t>((a >> shift) & 0b11);
    uint8_t bb = static_cast<uint8_t>((b >> shift) & 0b11);
    if (aa != bb) positions.push_back(i);
  }
  return positions;
}

const std::vector<SiteRecord> &candidate_sites(const GenomeIndex &index, const EngineParams &params) {
  if (params.search_backend == SearchBackend::BruteForce) {
    return index.sites();
  }
  throw std::runtime_error("SearchBackend::FMIndex is not implemented yet");
}

std::vector<OffTargetHit> run_cpu_engine(const GenomeIndex &index,
                                         const Guide &guide,
                                         const EngineParams &params) {
  auto wall_total_start = std::chrono::steady_clock::now();
  const auto &meta = index.meta();
  const auto &sites = candidate_sites(index, params);
  auto t1_start = std::chrono::steady_clock::now();
  EncodedGuide eg = encode_guide(guide, meta.guide_length);
  auto t1_end = std::chrono::steady_clock::now();

  std::vector<OffTargetHit> hits;
  hits.reserve(1024);

  auto t2_start = std::chrono::steady_clock::now();
  for (const auto &site : sites) {
    uint8_t mm = hamming_distance_2bit(eg.bits, site.seq_bits, meta.guide_length);
    if (mm > params.max_mismatches) continue;

    std::vector<uint8_t> mm_positions;
    if (params.score_params.model != ScoreModel::Hamming) {
      mm_positions = mismatch_positions(eg.bits, site.seq_bits, meta.guide_length);
    }

    float score = 0.0f;
    switch (params.score_params.model) {
      case ScoreModel::Hamming:
        score = score_mismatch_count(mm);
        break;
      case ScoreModel::MIT:
        score = score_mit(mm_positions, meta.guide_length);
        break;
      case ScoreModel::CFD:
        score = score_cfd_bits(eg.bits, site.seq_bits, meta.guide_length);
        break;
    }

    OffTargetHit hit;
    hit.guide_name = guide.name;
    hit.chrom_id = site.chrom_id;
    hit.pos = site.pos;
    hit.strand = site.strand == 0 ? '+' : '-';
    hit.mismatches = mm;
    hit.score = score;
    hits.push_back(hit);
  }
  auto t2_end = std::chrono::steady_clock::now();

  double s1 = std::chrono::duration<double>(t1_end - t1_start).count();
  double s2 = std::chrono::duration<double>(t2_end - t2_start).count();
  uint64_t candidates = static_cast<uint64_t>(sites.size());
  if (timing_enabled()) {
    log_stat("cpu.stage1", BenchStat{candidates, s1});
    log_stat("cpu.stage2", BenchStat{candidates, s2});
    double total = std::chrono::duration<double>(t2_end - wall_total_start).count();
    log_stat("cpu.total", BenchStat{candidates, total});
  }

  return hits;
}

#ifdef CRISPR_GPU_ENABLE_CUDA

} // namespace

void check_cuda(cudaError_t err, const char *msg) {
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string(msg) + ": " + cudaGetErrorString(err));
  }
}

namespace detail {

struct GpuContext {
  SiteRecord *d_sites{nullptr};
  DeviceHit *d_hits{nullptr};
  uint32_t *d_count{nullptr};
  uint32_t *h_count{nullptr};              // pinned
  DeviceHit *h_hits{nullptr};      // pinned
  size_t capacity_sites{0};
  size_t capacity_hits{0};
  const SiteRecord *host_sites_ptr{nullptr};
  size_t host_sites_size{0};
  cudaStream_t stream{nullptr};

  ~GpuContext() {
    if (d_sites) cudaFree(d_sites);
    if (d_hits) cudaFree(d_hits);
    if (d_count) cudaFree(d_count);
    if (h_count) cudaFreeHost(h_count);
    if (h_hits) cudaFreeHost(h_hits);
    if (stream) cudaStreamDestroy(stream);
  }

  void ensure_stream() {
    if (!stream) {
      check_cuda(cudaStreamCreate(&stream), "cudaStreamCreate");
    }
  }

  void ensure_capacity(size_t sites, size_t hits) {
    ensure_stream();
    if (sites > capacity_sites) {
      if (d_sites) cudaFree(d_sites);
      check_cuda(cudaMalloc(&d_sites, sites * sizeof(SiteRecord)), "cudaMalloc d_sites");
      capacity_sites = sites;
      // force re-upload on next call
      host_sites_ptr = nullptr;
      host_sites_size = 0;
    }
    if (hits > capacity_hits) {
      if (d_hits) cudaFree(d_hits);
      if (h_hits) cudaFreeHost(h_hits);
      check_cuda(cudaMalloc(&d_hits, hits * sizeof(DeviceHit)), "cudaMalloc d_hits");
      check_cuda(cudaHostAlloc(reinterpret_cast<void **>(&h_hits),
                               hits * sizeof(DeviceHit), cudaHostAllocDefault),
                 "cudaHostAlloc h_hits");
      capacity_hits = hits;
    }
    if (!d_count) {
      check_cuda(cudaMalloc(&d_count, sizeof(uint32_t)), "cudaMalloc d_count");
    }
    if (!h_count) {
      check_cuda(cudaHostAlloc(reinterpret_cast<void **>(&h_count), sizeof(uint32_t),
                               cudaHostAllocDefault),
                 "cudaHostAlloc h_count");
    }
  }

  void ensure_sites_uploaded(const std::vector<SiteRecord> &sites) {
    if (host_sites_ptr == sites.data() && host_sites_size == sites.size()) {
      return; // already resident
    }
    ensure_capacity(sites.size(), sites.size());
    check_cuda(cudaMemcpyAsync(d_sites, sites.data(),
                               sites.size() * sizeof(SiteRecord),
                               cudaMemcpyHostToDevice, stream),
               "cudaMemcpyAsync sites");
    host_sites_ptr = sites.data();
    host_sites_size = sites.size();
    check_cuda(cudaStreamSynchronize(stream), "cudaStreamSync after sites upload");
  }
};

struct EngineGpuState {
  detail::GpuContext ctx;
  bool sites_ready{false};
};

} // namespace detail

namespace {

using detail::EngineGpuState;
using detail::GpuContext;

void cuda_warmup_internal() {
  static detail::GpuContext ctx; // isolated warmup context
  const uint64_t dummy_bits = 0;
  // allocate minimal buffers once
  ctx.ensure_capacity(1, 1);
  SiteRecord dummy{};
  dummy.seq_bits = dummy_bits;
  dummy.chrom_id = 0;
  dummy.pos = 0;
  dummy.strand = 0;
  std::vector<SiteRecord> host_dummy(1, dummy);
  ctx.ensure_sites_uploaded(host_dummy);
  check_cuda(cudaMemsetAsync(ctx.d_count, 0, sizeof(uint32_t), ctx.stream), "warmup memset");
  detail::launch_off_target_kernel(ctx.d_sites, 1, dummy_bits, 0, 20,
                                   ScoreParams{}, ctx.d_hits, ctx.d_count, ctx.stream);
  check_cuda(cudaStreamSynchronize(ctx.stream), "warmup sync");
}

std::vector<OffTargetHit> run_gpu_engine(const GenomeIndex &index,
                                         const Guide &guide,
                                         const EngineParams &params,
                                         EngineGpuState *state);

std::vector<OffTargetHit> run_gpu_engine_batch(const GenomeIndex &index,
                                               const std::vector<Guide> &guides,
                                               const EngineParams &params,
                                               EngineGpuState *state) {
  if (guides.empty()) return {};

  const auto &sites = candidate_sites(index, params);
  uint32_t num_sites = static_cast<uint32_t>(sites.size());
  if (num_sites == 0) return {};

  static detail::GpuContext fallback;
  detail::GpuContext &ctx = state ? state->ctx : fallback;
  ctx.ensure_capacity(num_sites, num_sites);
  if (!state || !state->sites_ready) {
    ScopedTimer t_upload("gpu.upload_sites", timing_enabled());
    ctx.ensure_sites_uploaded(sites);
    if (state) state->sites_ready = true;
  }

  auto wall_total_start = std::chrono::steady_clock::now();
  auto t1_start = std::chrono::steady_clock::now();
  std::vector<EncodedGuide> encoded;
  encoded.reserve(guides.size());
  for (const auto &g : guides) {
    encoded.push_back(encode_guide(g, index.meta().guide_length));
  }
  auto t1_end = std::chrono::steady_clock::now();

  double stage2_accum = 0.0;

  std::vector<OffTargetHit> all_hits;
  for (size_t gi = 0; gi < guides.size(); ++gi) {
    // zero count per guide
    {
      ScopedTimer t_zero("gpu.zero", timing_enabled());
      check_cuda(cudaMemsetAsync(ctx.d_count, 0, sizeof(uint32_t), ctx.stream), "cudaMemsetAsync count");
    }

    auto t2_start = std::chrono::steady_clock::now();
    {
      ScopedTimer t_kernel("gpu.kernel", timing_enabled());
      detail::launch_off_target_kernel(ctx.d_sites, num_sites, encoded[gi].bits, params.max_mismatches,
                                       index.meta().guide_length, params.score_params,
                                       ctx.d_hits, ctx.d_count, ctx.stream);
    }

    {
      ScopedTimer t_count("gpu.copy_count", timing_enabled());
      check_cuda(cudaMemcpyAsync(ctx.h_count, ctx.d_count, sizeof(uint32_t),
                                 cudaMemcpyDeviceToHost, ctx.stream),
                 "cudaMemcpyAsync count back");
    }
    check_cuda(cudaStreamSynchronize(ctx.stream), "cudaStreamSynchronize after kernel");

    uint32_t host_count = *ctx.h_count;
    host_count = std::min<uint32_t>(host_count, num_sites);

    if (host_count > 0) {
      ScopedTimer t_hits("gpu.copy_hits", timing_enabled());
      check_cuda(cudaMemcpyAsync(ctx.h_hits, ctx.d_hits,
                                 host_count * sizeof(detail::DeviceHit),
                                 cudaMemcpyDeviceToHost, ctx.stream),
                 "cudaMemcpyAsync hits back");
      check_cuda(cudaStreamSynchronize(ctx.stream), "cudaStreamSynchronize after hits");
    }
    auto t2_end = std::chrono::steady_clock::now();
    stage2_accum += std::chrono::duration<double>(t2_end - t2_start).count();

    all_hits.reserve(all_hits.size() + host_count);
    for (uint32_t i = 0; i < host_count; ++i) {
      const auto &h = ctx.h_hits[i];
      const auto &site = sites[h.site_index];
      OffTargetHit o;
      o.guide_name = guides[gi].name;
      o.chrom_id = site.chrom_id;
      o.pos = site.pos;
      o.strand = site.strand == 0 ? '+' : '-';
      o.mismatches = h.mismatches;
      o.score = h.score;
      all_hits.push_back(o);
    }
  }

  if (timing_enabled()) {
    double s1 = std::chrono::duration<double>(t1_end - t1_start).count();
    double s2 = stage2_accum;
    uint64_t candidates = static_cast<uint64_t>(sites.size()) * static_cast<uint64_t>(guides.size());
    double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_total_start).count();
    log_stat("gpu.stage1", BenchStat{candidates, s1});
    log_stat("gpu.stage2", BenchStat{candidates, s2});
    log_stat("gpu.total", BenchStat{candidates, total});
  }

  return all_hits;
}
std::vector<OffTargetHit> run_gpu_engine(const GenomeIndex &index,
                                         const Guide &guide,
                                         const EngineParams &params,
                                         EngineGpuState *state) {
  auto wall_total_start = std::chrono::steady_clock::now();
  const auto &sites = candidate_sites(index, params);
  uint32_t num_sites = static_cast<uint32_t>(sites.size());
  if (num_sites == 0) return {};

  auto t1_start = std::chrono::steady_clock::now();
  EncodedGuide eg = encode_guide(guide, index.meta().guide_length);
  auto t1_end = std::chrono::steady_clock::now();

  static detail::GpuContext fallback;
  detail::GpuContext &ctx = state ? state->ctx : fallback;
  ctx.ensure_capacity(num_sites, num_sites);
  if (!state || !state->sites_ready) {
    ScopedTimer t_upload("gpu.upload_sites", timing_enabled());
    ctx.ensure_sites_uploaded(sites);
    if (state) state->sites_ready = true;
  }

  // zero count
  {
    ScopedTimer t_zero("gpu.zero", timing_enabled());
    check_cuda(cudaMemsetAsync(ctx.d_count, 0, sizeof(uint32_t), ctx.stream), "cudaMemsetAsync count");
  }

  auto t2_start = std::chrono::steady_clock::now();
  {
    ScopedTimer t_kernel("gpu.kernel", timing_enabled());
    detail::launch_off_target_kernel(ctx.d_sites, num_sites, eg.bits, params.max_mismatches,
                                     index.meta().guide_length, params.score_params,
                                     ctx.d_hits, ctx.d_count, ctx.stream);
  }

  {
    ScopedTimer t_count("gpu.copy_count", timing_enabled());
    check_cuda(cudaMemcpyAsync(ctx.h_count, ctx.d_count, sizeof(uint32_t),
                               cudaMemcpyDeviceToHost, ctx.stream),
               "cudaMemcpyAsync count back");
  }
  check_cuda(cudaStreamSynchronize(ctx.stream), "cudaStreamSynchronize after kernel");

  uint32_t host_count = *ctx.h_count;
  host_count = std::min<uint32_t>(host_count, num_sites);

  if (host_count > 0) {
    ScopedTimer t_hits("gpu.copy_hits", timing_enabled());
    check_cuda(cudaMemcpyAsync(ctx.h_hits, ctx.d_hits,
                               host_count * sizeof(detail::DeviceHit),
                               cudaMemcpyDeviceToHost, ctx.stream),
               "cudaMemcpyAsync hits back");
    check_cuda(cudaStreamSynchronize(ctx.stream), "cudaStreamSynchronize after hits");
  }
  auto t2_end = std::chrono::steady_clock::now();

  std::vector<OffTargetHit> out;
  out.reserve(host_count);
  for (uint32_t i = 0; i < host_count; ++i) {
    const auto &h = ctx.h_hits[i];
    const auto &site = sites[h.site_index];
    OffTargetHit o;
    o.guide_name = guide.name;
    o.chrom_id = site.chrom_id;
    o.pos = site.pos;
    o.strand = site.strand == 0 ? '+' : '-';
    o.mismatches = h.mismatches;
    o.score = h.score;
    out.push_back(o);
  }

  if (timing_enabled()) {
    double s1 = std::chrono::duration<double>(t1_end - t1_start).count();
    double s2 = std::chrono::duration<double>(t2_end - t2_start).count();
    uint64_t candidates = static_cast<uint64_t>(sites.size());
    log_stat("gpu.stage1", BenchStat{candidates, s1});
    log_stat("gpu.stage2", BenchStat{candidates, s2});
    double total = std::chrono::duration<double>(t2_end - wall_total_start).count();
    log_stat("gpu.total", BenchStat{candidates, total});
  }

  return out;
}
#endif // CRISPR_GPU_ENABLE_CUDA

} // namespace

OffTargetEngine::OffTargetEngine(const GenomeIndex &index, EngineParams params)
    : index_(index), params_(params) {
  if (params_.max_mismatches > index.meta().guide_length) {
    throw std::runtime_error("max_mismatches cannot exceed guide length");
  }
#ifndef CRISPR_GPU_ENABLE_CUDA
  params_.backend = Backend::CPU;
#endif
  const char *env_backend = std::getenv("CRISPR_GPU_BACKEND");
  if (env_backend) {
    std::string v(env_backend);
    std::transform(v.begin(), v.end(), v.begin(), ::tolower);
    if (v == "cpu") params_.backend = Backend::CPU;
    if (v == "gpu") params_.backend = Backend::GPU;
  }
#ifdef CRISPR_GPU_ENABLE_CUDA
  int dev_count = 0;
  if (params_.backend == Backend::GPU) {
    cudaError_t err = cudaGetDeviceCount(&dev_count);
    if (err != cudaSuccess || dev_count == 0) {
      params_.backend = Backend::CPU;
    }
  }
  if (params_.backend == Backend::GPU) {
    gpu_state_ = std::make_unique<EngineGpuState>();
  }
#endif
}

OffTargetEngine::~OffTargetEngine() = default;

std::vector<OffTargetHit> OffTargetEngine::score_guide(const Guide &guide) const {
  ensure_default_tables_loaded();
  if (params_.backend == Backend::GPU) {
#ifdef CRISPR_GPU_ENABLE_CUDA
    if (!gpu_state_) gpu_state_ = std::make_unique<EngineGpuState>();
    return run_gpu_engine(index_, guide, params_, gpu_state_.get());
#else
    return run_cpu_engine(index_, guide, params_);
#endif
  }
  return run_cpu_engine(index_, guide, params_);
}

std::vector<OffTargetHit> OffTargetEngine::score_guides(const std::vector<Guide> &guides) const {
  auto t_start = std::chrono::steady_clock::now();

  std::vector<OffTargetHit> all;
  ensure_default_tables_loaded();
#ifdef CRISPR_GPU_ENABLE_CUDA
  if (params_.backend == Backend::GPU) {
    if (!gpu_state_) gpu_state_ = std::make_unique<EngineGpuState>();
    all = run_gpu_engine_batch(index_, guides, params_, gpu_state_.get());
  } else {
    for (const auto &g : guides) {
      auto h = run_cpu_engine(index_, g, params_);
      all.insert(all.end(), h.begin(), h.end());
    }
  }
#else
  for (const auto &g : guides) {
    auto h = run_cpu_engine(index_, g, params_);
    all.insert(all.end(), h.begin(), h.end());
  }
#endif

  if (timing_enabled()) {
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    uint64_t candidates = static_cast<uint64_t>(index_.sites().size()) * static_cast<uint64_t>(guides.size());
    const char *label = params_.backend == Backend::GPU ? "gpu.guides" : "cpu.guides";
    log_stat(label, BenchStat{candidates, seconds});
  }

  return all;
}

bool cuda_available() {
#ifdef CRISPR_GPU_ENABLE_CUDA
  int dev_count = 0;
  if (cudaGetDeviceCount(&dev_count) != cudaSuccess) return false;
  return dev_count > 0;
#else
  return false;
#endif
}

void cuda_warmup() {
#ifdef CRISPR_GPU_ENABLE_CUDA
  try {
    cuda_available(); // triggers context creation
    cuda_warmup_internal(); // call the internal helper
  } catch (...) {
  }
#endif
}

void log_stat(const char *label, const BenchStat &s) {
  std::fprintf(stderr, "[stat] %s: candidates=%llu time=%.6f cgct=%.3e\n",
               label,
               static_cast<unsigned long long>(s.candidates),
               s.seconds,
               s.cgct());
}

} // namespace crispr_gpu
