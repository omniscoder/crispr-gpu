#include <catch2/catch_all.hpp>
#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/engine.hpp"
#include <tuple>
#include <algorithm>
#include <fstream>
#include <unistd.h>
#ifdef CRISPR_GPU_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

using namespace crispr_gpu;

#ifdef CRISPR_GPU_ENABLE_CUDA

TEST_CASE("GPU and CPU engines agree on toy genome", "[gpu]") {
  IndexParams params;
  params.guide_length = 4;
  params.pam = "NGG";
  params.both_strands = true;

#ifdef CRISPR_GPU_ENABLE_CUDA
  int device_count = 0;
  cudaError_t err = cudaGetDeviceCount(&device_count);
  if (err != cudaSuccess || device_count == 0) {
    SUCCEED("No CUDA device available");
    return;
  }
#endif

  auto idx = GenomeIndex::build(std::string(TEST_DATA_DIR) + "/toy.fa", params);

  Guide g;
  g.name = "g1";
  g.sequence = "AAAA";
  g.pam = "NGG";

  EngineParams cpu_params;
  cpu_params.max_mismatches = 4;
  cpu_params.backend = Backend::CPU;

  EngineParams gpu_params = cpu_params;
  gpu_params.backend = Backend::GPU;

  OffTargetEngine cpu_engine(idx, cpu_params);
  OffTargetEngine gpu_engine(idx, gpu_params);

  auto cpu_hits = cpu_engine.score_guide(g);
  auto gpu_hits = gpu_engine.score_guide(g);

  REQUIRE(cpu_hits.size() == gpu_hits.size());
  // Sort both vectors to deterministic order
  auto key = [](const OffTargetHit &h) {
    return std::tuple<uint32_t, uint32_t, char, uint8_t>(h.chrom_id, h.pos, h.strand, h.mismatches);
  };
  std::sort(cpu_hits.begin(), cpu_hits.end(), [&](const auto &a, const auto &b){ return key(a) < key(b);});
  std::sort(gpu_hits.begin(), gpu_hits.end(), [&](const auto &a, const auto &b){ return key(a) < key(b);});
  for (size_t i = 0; i < cpu_hits.size(); ++i) {
    REQUIRE(cpu_hits[i].chrom_id == gpu_hits[i].chrom_id);
    REQUIRE(cpu_hits[i].pos == gpu_hits[i].pos);
    REQUIRE(cpu_hits[i].strand == gpu_hits[i].strand);
    REQUIRE(cpu_hits[i].mismatches == gpu_hits[i].mismatches);
    REQUIRE(cpu_hits[i].score == Catch::Approx(gpu_hits[i].score).epsilon(1e-6));
  }
}

#ifdef CRISPR_GPU_ENABLE_CUDA

TEST_CASE("GPU and CPU agree on MIT/CFD scoring placeholders", "[gpu]") {
  int device_count = 0;
  cudaError_t err = cudaGetDeviceCount(&device_count);
  if (err != cudaSuccess || device_count == 0) {
    SUCCEED("No CUDA device available");
    return;
  }

  char fname[] = "/tmp/crispr_gpu_mitXXXXXX";
  int fd = mkstemp(fname);
  REQUIRE(fd != -1);
  close(fd);
  std::string fasta = fname;
  {
    std::ofstream out(fasta);
    out << ">chr1\n";
    out << "AAAA"     // perfect
        << "ACAA"     // 1 mismatch (pos1)
        << "AACA";    // 1 mismatch (pos2)
    out << "\n";
  }

  IndexParams params;
  params.guide_length = 4;
  params.pam = "NNN"; // accept everywhere to keep builder simple
  params.both_strands = false;
  auto idx = GenomeIndex::build(fasta, params);

  Guide g{"g1", "AAAA", "NNN"};

  EngineParams cpu_params;
  cpu_params.max_mismatches = 4;
  cpu_params.backend = Backend::CPU;
  cpu_params.score_params.model = ScoreModel::MIT;

  EngineParams gpu_params = cpu_params;
  gpu_params.backend = Backend::GPU;

  OffTargetEngine cpu_engine(idx, cpu_params);
  OffTargetEngine gpu_engine(idx, gpu_params);

  auto cpu_hits = cpu_engine.score_guide(g);
  auto gpu_hits = gpu_engine.score_guide(g);

  REQUIRE(cpu_hits.size() == gpu_hits.size());
  auto key = [](const OffTargetHit &h) {
    return std::tuple<uint32_t, uint32_t, char, uint8_t>(h.chrom_id, h.pos, h.strand, h.mismatches);
  };
  std::sort(cpu_hits.begin(), cpu_hits.end(), [&](const auto &a, const auto &b){ return key(a) < key(b);});
  std::sort(gpu_hits.begin(), gpu_hits.end(), [&](const auto &a, const auto &b){ return key(a) < key(b);});
  for (size_t i = 0; i < cpu_hits.size(); ++i) {
    REQUIRE(cpu_hits[i].chrom_id == gpu_hits[i].chrom_id);
    REQUIRE(cpu_hits[i].pos == gpu_hits[i].pos);
    REQUIRE(cpu_hits[i].strand == gpu_hits[i].strand);
    REQUIRE(cpu_hits[i].mismatches == gpu_hits[i].mismatches);
    REQUIRE(cpu_hits[i].score == Catch::Approx(gpu_hits[i].score).epsilon(1e-6));
  }
}

#endif // CRISPR_GPU_ENABLE_CUDA

#endif // CRISPR_GPU_ENABLE_CUDA
