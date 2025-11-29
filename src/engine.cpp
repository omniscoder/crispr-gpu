#include "crispr_gpu/engine.hpp"
#include "crispr_gpu/types.hpp"
#include "crispr_gpu/scoring.hpp"

#include <stdexcept>
#include <algorithm>
#include <numeric>

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

std::vector<OffTargetHit> run_cpu_engine(const GenomeIndex &index,
                                         const Guide &guide,
                                         const EngineParams &params) {
  const auto &meta = index.meta();
  EncodedGuide eg = encode_guide(guide, meta.guide_length);

  std::vector<OffTargetHit> hits;
  hits.reserve(1024);

  for (const auto &site : index.sites()) {
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

  return hits;
}

#ifdef CRISPR_GPU_ENABLE_CUDA

void check_cuda(cudaError_t err, const char *msg) {
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string(msg) + ": " + cudaGetErrorString(err));
  }
}

std::vector<OffTargetHit> run_gpu_engine(const GenomeIndex &index,
                                         const Guide &guide,
                                         const EngineParams &params) {
  const auto &sites = index.sites();
  uint32_t num_sites = static_cast<uint32_t>(sites.size());
  if (num_sites == 0) return {};

  EncodedGuide eg = encode_guide(guide, index.meta().guide_length);

  SiteRecord *d_sites = nullptr;
  detail::DeviceHit *d_hits = nullptr;
  uint32_t *d_count = nullptr;
  uint32_t max_hits = num_sites;

  check_cuda(cudaMalloc(&d_sites, num_sites * sizeof(SiteRecord)), "cudaMalloc d_sites");
  check_cuda(cudaMalloc(&d_hits, max_hits * sizeof(detail::DeviceHit)), "cudaMalloc d_hits");
  check_cuda(cudaMalloc(&d_count, sizeof(uint32_t)), "cudaMalloc d_count");

  check_cuda(cudaMemcpy(d_sites, sites.data(), num_sites * sizeof(SiteRecord), cudaMemcpyHostToDevice),
             "cudaMemcpy sites");
  uint32_t zero = 0;
  check_cuda(cudaMemcpy(d_count, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice), "cudaMemcpy count");

  detail::launch_off_target_kernel(d_sites, num_sites, eg.bits, params.max_mismatches,
                                   index.meta().guide_length, params.score_params, d_hits, d_count, 0);
  check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

  uint32_t host_count = 0;
  check_cuda(cudaMemcpy(&host_count, d_count, sizeof(uint32_t), cudaMemcpyDeviceToHost),
             "cudaMemcpy count back");
  host_count = std::min(host_count, max_hits);

  std::vector<detail::DeviceHit> device_hits(host_count);
  if (host_count > 0) {
    check_cuda(cudaMemcpy(device_hits.data(), d_hits, host_count * sizeof(detail::DeviceHit),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy hits back");
  }

  cudaFree(d_sites);
  cudaFree(d_hits);
  cudaFree(d_count);

  std::vector<OffTargetHit> out;
  out.reserve(host_count);
  for (const auto &h : device_hits) {
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
  cudaError_t err = cudaGetDeviceCount(&dev_count);
  if (params_.backend == Backend::GPU && (err != cudaSuccess || dev_count == 0)) {
    params_.backend = Backend::CPU;
  }
#endif
}

std::vector<OffTargetHit> OffTargetEngine::score_guide(const Guide &guide) const {
  ensure_default_tables_loaded();
  if (params_.backend == Backend::GPU) {
#ifdef CRISPR_GPU_ENABLE_CUDA
    return run_gpu_engine(index_, guide, params_);
#else
    return run_cpu_engine(index_, guide, params_);
#endif
  }
  return run_cpu_engine(index_, guide, params_);
}

std::vector<OffTargetHit> OffTargetEngine::score_guides(const std::vector<Guide> &guides) const {
  std::vector<OffTargetHit> all;
  for (const auto &g : guides) {
    auto h = score_guide(g);
    all.insert(all.end(), h.begin(), h.end());
  }
  return all;
}

} // namespace crispr_gpu
