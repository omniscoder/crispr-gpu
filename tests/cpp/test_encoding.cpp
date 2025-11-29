#include <catch2/catch_all.hpp>
#include "crispr_gpu/types.hpp"

using namespace crispr_gpu;

TEST_CASE("2bit encoding round-trip") {
  std::string seq = "ACGTACGTACGTACGTACGT"; // 20 nt
  auto bits = encode_sequence_2bit(seq);
  auto decoded = decode_sequence_2bit(bits, static_cast<uint8_t>(seq.size()));
  REQUIRE(decoded == seq);
}

TEST_CASE("Hamming distance") {
  std::string a = "AAAAAAAAAAAAAAAAAAAA";
  std::string b = "AAAACAAAAAAAAAAAAAAA"; // 1 mismatch
  auto bits_a = encode_sequence_2bit(a);
  auto bits_b = encode_sequence_2bit(b);
  auto dist = hamming_distance_2bit(bits_a, bits_b, 20);
  REQUIRE(dist == 1);
}
