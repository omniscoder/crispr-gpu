#include "crispr_gpu/types.hpp"

#include <algorithm>
#include <stdexcept>

namespace crispr_gpu {

uint8_t base_to_bits(char b) {
  switch (b) {
    case 'A': case 'a': return 0b00;
    case 'C': case 'c': return 0b01;
    case 'G': case 'g': return 0b10;
    case 'T': case 't': return 0b11;
    default:
      throw std::runtime_error(std::string("Invalid base: ") + b);
  }
}

char bits_to_base(uint8_t bits) {
  bits &= 0b11;
  switch (bits) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    default: return 'N';
  }
}

uint64_t encode_sequence_2bit(const std::string &seq) {
  if (seq.size() > 32) {
    throw std::runtime_error("encode_sequence_2bit supports sequences up to 32 bases");
  }
  uint64_t bits = 0;
  for (char c : seq) {
    bits = (bits << 2) | base_to_bits(c);
  }
  return bits;
}

std::string decode_sequence_2bit(uint64_t bits, uint8_t length) {
  std::string out(length, 'A');
  for (int i = length - 1; i >= 0; --i) {
    out[static_cast<size_t>(i)] = bits_to_base(static_cast<uint8_t>(bits & 0b11));
    bits >>= 2;
  }
  return out;
}

std::string revcomp(const std::string &seq) {
  std::string out(seq.size(), 'N');
  for (size_t i = 0; i < seq.size(); ++i) {
    char c = seq[seq.size() - 1 - i];
    switch (c) {
      case 'A': case 'a': out[i] = 'T'; break;
      case 'C': case 'c': out[i] = 'G'; break;
      case 'G': case 'g': out[i] = 'C'; break;
      case 'T': case 't': out[i] = 'A'; break;
      default: out[i] = 'N'; break;
    }
  }
  return out;
}

uint8_t hamming_distance_2bit(uint64_t a, uint64_t b, uint8_t length) {
  uint8_t used_bits = static_cast<uint8_t>(length * 2);
  uint64_t mask = (used_bits == 64) ? ~0ULL : ((1ULL << used_bits) - 1);
  uint64_t x = (a ^ b) & mask;
  uint64_t mism = (x | (x >> 1)) & 0x5555555555555555ULL;
  return static_cast<uint8_t>(__builtin_popcountll(mism));
}

} // namespace crispr_gpu
