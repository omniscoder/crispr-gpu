#include <catch2/catch_all.hpp>
#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/engine.hpp"

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

using namespace crispr_gpu;

TEST_CASE("Index build and score on toy genome") {
  char fname[] = "/tmp/crispr_gpu_toyXXXXXX";
  int fd = mkstemp(fname);
  REQUIRE(fd != -1);
  close(fd);
  std::string fasta = fname;
  std::ofstream out(fasta);
  out << ">chr1\n";
  out << "AAAAGGGAAAA\n"; // protospacer AAAA, PAM GGG
  out.close();

  IndexParams params;
  params.guide_length = 4;
  params.pam = "NGG";
  params.both_strands = true;

  auto idx = GenomeIndex::build(fasta, params);
  REQUIRE(idx.sites().size() >= 1);
  REQUIRE(idx.meta().guide_length == 4);

  Guide g{"g1", "AAAA", "NGG"};
  OffTargetEngine engine(idx, {});
  auto hits = engine.score_guide(g);
  REQUIRE(!hits.empty());
  auto hit0 = hits.front();
  REQUIRE(hit0.mismatches == 0);
}
