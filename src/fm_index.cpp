#include "crispr_gpu/fm_index.hpp"

#include <fstream>
#include <stdexcept>
#include <algorithm>

namespace crispr_gpu {

namespace {

inline uint8_t dna_to_sym(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 0; // treat unknown as A for now
  }
}

} // namespace

FmIndex build_fm_index(const std::vector<std::string> &protospacers,
                       const std::vector<FmIndex::Locus> &loci,
                       const FmIndexHeader &header) {
  if (protospacers.size() != loci.size()) {
    throw std::runtime_error("protospacers and loci size mismatch");
  }
  // Placeholder: real FM build to be implemented.
  FmIndex idx;
  idx.header = header;
  idx.loci = loci;
  idx.header.text_length = static_cast<uint64_t>(protospacers.size()) * header.guide_length;
  return idx;
}

void save_fm_index(const FmIndex &index, const std::string &path) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot open fm index for write: " + path);
  // Very simple placeholder serialization; to be upgraded to a stable layout.
  out.write(reinterpret_cast<const char *>(&index.header.version), sizeof(index.header.version));
  out.write(reinterpret_cast<const char *>(&index.header.guide_length), sizeof(index.header.guide_length));
  uint64_t pam_len = index.header.pam.size();
  out.write(reinterpret_cast<const char *>(&pam_len), sizeof(pam_len));
  out.write(index.header.pam.data(), pam_len);
  out.write(reinterpret_cast<const char *>(&index.header.both_strands), sizeof(index.header.both_strands));
  out.write(reinterpret_cast<const char *>(&index.header.text_length), sizeof(index.header.text_length));
  out.write(reinterpret_cast<const char *>(&index.header.sampled_sa_interval), sizeof(index.header.sampled_sa_interval));

  auto write_vec = [&out](const auto &v) {
    uint64_t n = v.size();
    out.write(reinterpret_cast<const char *>(&n), sizeof(n));
    out.write(reinterpret_cast<const char *>(v.data()), n * sizeof(typename std::decay_t<decltype(v)>::value_type));
  };

  write_vec(index.occ_A);
  write_vec(index.occ_C);
  write_vec(index.occ_G);
  write_vec(index.occ_T);
  write_vec(index.sa_samples);
  write_vec(index.loci);
}

FmIndex load_fm_index(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open fm index for read: " + path);
  FmIndex idx;
  in.read(reinterpret_cast<char *>(&idx.header.version), sizeof(idx.header.version));
  in.read(reinterpret_cast<char *>(&idx.header.guide_length), sizeof(idx.header.guide_length));
  uint64_t pam_len = 0;
  in.read(reinterpret_cast<char *>(&pam_len), sizeof(pam_len));
  idx.header.pam.resize(pam_len, '\0');
  in.read(&idx.header.pam[0], pam_len);
  in.read(reinterpret_cast<char *>(&idx.header.both_strands), sizeof(idx.header.both_strands));
  in.read(reinterpret_cast<char *>(&idx.header.text_length), sizeof(idx.header.text_length));
  in.read(reinterpret_cast<char *>(&idx.header.sampled_sa_interval), sizeof(idx.header.sampled_sa_interval));

  auto read_vec = [&in](auto &v) {
    uint64_t n = 0; in.read(reinterpret_cast<char *>(&n), sizeof(n));
    using T = typename std::decay_t<decltype(v)>::value_type;
    v.resize(n);
    in.read(reinterpret_cast<char *>(v.data()), n * sizeof(T));
  };

  read_vec(idx.occ_A);
  read_vec(idx.occ_C);
  read_vec(idx.occ_G);
  read_vec(idx.occ_T);
  read_vec(idx.sa_samples);
  read_vec(idx.loci);
  return idx;
}

std::vector<uint32_t> fm_search_hamming(const FmIndex &index,
                                        const std::string &pattern,
                                        uint8_t max_mismatches) {
  (void)index; (void)pattern; (void)max_mismatches;
  // Placeholder: return empty until implemented.
  return {};
}

} // namespace crispr_gpu
