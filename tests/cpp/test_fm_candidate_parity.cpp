#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/engine.hpp"

#include <catch2/catch.hpp>
#include <algorithm>

using namespace crispr_gpu;

namespace {

std::vector<Guide> load_guides(const std::string &path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open guides file");
  std::vector<Guide> guides;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    Guide g; ss >> g.name >> g.sequence >> g.pam;
    guides.push_back(g);
  }
  return guides;
}

auto sort_sites(std::vector<SiteRecord> v) {
  std::sort(v.begin(), v.end(), [](const SiteRecord &a, const SiteRecord &b) {
    if (a.chrom_id != b.chrom_id) return a.chrom_id < b.chrom_id;
    if (a.pos != b.pos) return a.pos < b.pos;
    return a.strand < b.strand;
  });
  return v;
}

} // namespace

TEST_CASE("FM backend matches brute for K=0 and K=1", "[fm][parity]") {
  const std::string fasta = "tests/data/fm_tiny.fa";
  const std::string guides_path = "tests/data/fm_tiny_guides.tsv";

  IndexParams ip; ip.guide_length = 12; ip.pam = "NGG"; ip.both_strands = true;
  auto idx = GenomeIndex::build(fasta, ip);
  auto guides = load_guides(guides_path);

  for (uint8_t K : {uint8_t(0), uint8_t(1)}) {
    EngineParams brute_params; brute_params.max_mismatches = K; brute_params.search_backend = SearchBackend::BruteForce; brute_params.backend = Backend::CPU;
    EngineParams fm_params = brute_params; fm_params.search_backend = SearchBackend::FMIndex;

    OffTargetEngine eng_brute(idx, brute_params);
    OffTargetEngine eng_fm(idx, fm_params);

    for (const auto &g : guides) {
      auto sb = eng_brute.score_guide(g);
      auto sf = eng_fm.score_guide(g);
      auto key_hits = [](const std::vector<OffTargetHit> &hits) {
        auto sorted = hits;
        std::sort(sorted.begin(), sorted.end(), [](const OffTargetHit &a, const OffTargetHit &b) {
          if (a.chrom_id != b.chrom_id) return a.chrom_id < b.chrom_id;
          if (a.pos != b.pos) return a.pos < b.pos;
          if (a.strand != b.strand) return a.strand < b.strand;
          if (a.mismatches != b.mismatches) return a.mismatches < b.mismatches;
          return a.guide_name < b.guide_name;
        });
        return sorted;
      };
      REQUIRE(key_hits(sb) == key_hits(sf));
    }
  }
}

