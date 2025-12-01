#include "crispr_gpu/fm_index.hpp"

#include <catch2/catch.hpp>

using namespace crispr_gpu;

TEST_CASE("FM exact search matches naive find", "[fm]") {
  std::string contig = "ACGTACGTACGT"; // length 12
  std::vector<std::string> protos = {contig};
  std::vector<FmIndex::Locus> loci = {{0, 0, 0}};
  FmIndexHeader hdr;
  hdr.guide_length = static_cast<uint8_t>(contig.size());
  hdr.occ_block = 4;
  hdr.sa_sample_rate = 2;

  auto fm = build_fm_index(protos, loci, hdr);

  std::string pattern = "CGT";
  auto rows = fm_search_hamming(fm, pattern, 0);
  std::vector<uint32_t> poses;
  for (auto r : rows) poses.push_back(fm_row_to_pos(fm, r));
  std::sort(poses.begin(), poses.end());

  std::vector<uint32_t> expected;
  for (size_t p = contig.find(pattern); p != std::string::npos; p = contig.find(pattern, p + 1)) {
    expected.push_back(static_cast<uint32_t>(p));
  }
  std::sort(expected.begin(), expected.end());
  REQUIRE(poses == expected);
}

TEST_CASE("FM mismatch search finds variants", "[fm]") {
  std::string contig = "ACGTACGTACGT";
  std::vector<std::string> protos = {contig};
  std::vector<FmIndex::Locus> loci = {{0, 0, 0}};
  FmIndexHeader hdr;
  hdr.guide_length = static_cast<uint8_t>(contig.size());
  hdr.occ_block = 4;
  hdr.sa_sample_rate = 2;
  auto fm = build_fm_index(protos, loci, hdr);

  std::string pattern = "ACGT"; // occurs at 0,4,8
  // Introduce one mismatch: AGGT should match positions with 1 mismatch.
  auto rows = fm_search_hamming(fm, "AGGT", 1);
  std::vector<uint32_t> poses;
  for (auto r : rows) poses.push_back(fm_row_to_pos(fm, r));
  std::sort(poses.begin(), poses.end());
  std::vector<uint32_t> expected = {0,4,8};
  REQUIRE(poses == expected);
}

