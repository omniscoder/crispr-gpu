#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include "crispr_gpu/types.hpp"

namespace crispr_gpu {

struct IndexParams {
  uint8_t guide_length{20};
  std::string pam{"NGG"};
  bool both_strands{true};
};

class GenomeIndex {
public:
  GenomeIndex() = default;

  static GenomeIndex build(const std::string &fasta_path, const IndexParams &params);
  static GenomeIndex load(const std::string &index_path);

  void save(const std::string &index_path) const;

  const std::vector<SiteRecord> &sites() const { return sites_; }
  const std::vector<ChromInfo> &chromosomes() const { return chroms_; }
  const IndexMeta &meta() const { return meta_; }

private:
  IndexMeta meta_{};
  std::vector<ChromInfo> chroms_;
  std::vector<SiteRecord> sites_;
};

} // namespace crispr_gpu
