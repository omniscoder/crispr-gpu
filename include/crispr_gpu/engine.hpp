#pragma once

#include <vector>
#include "crispr_gpu/types.hpp"
#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/scoring.hpp"

namespace crispr_gpu {

namespace detail {
struct DeviceHit {
  uint32_t site_index;
  uint8_t mismatches;
  float score;
};
} // namespace detail

class OffTargetEngine {
public:
  OffTargetEngine(const GenomeIndex &index, EngineParams params = {});

  std::vector<OffTargetHit> score_guide(const Guide &guide) const;
  std::vector<OffTargetHit> score_guides(const std::vector<Guide> &guides) const;

private:
  const GenomeIndex &index_;
  EngineParams params_;
};

// True if this build has CUDA enabled and at least one device is present.
bool cuda_available();

} // namespace crispr_gpu
