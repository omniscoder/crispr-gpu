#include "crispr_gpu/fm_index.hpp"

#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdexcept>

namespace crispr_gpu {

namespace {

inline uint8_t dna_to_sym(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4; // N / sentinel
  }
}

// Simple O(n log n) SA for modest texts (good enough for first pass).
std::vector<uint32_t> build_sa(const std::vector<uint8_t> &text) {
  uint32_t n = static_cast<uint32_t>(text.size());
  std::vector<uint32_t> sa(n);
  std::iota(sa.begin(), sa.end(), 0u);
  std::sort(sa.begin(), sa.end(), [&](uint32_t a, uint32_t b) {
    return std::lexicographical_compare(text.begin() + a, text.end(), text.begin() + b, text.end());
  });
  return sa;
}

uint32_t occ(const FmIndex &fm, uint8_t c, uint32_t i) {
  if (i == 0) return 0;
  const uint32_t block = fm.header.occ_block;
  uint32_t block_idx = i / block;
  uint32_t base = fm.occ_checkpoints[block_idx][c];
  uint32_t start = block_idx * block;
  uint32_t end = i;
  for (uint32_t k = start; k < end && k < fm.text_len; ++k) {
    if (fm.bwt[k] == c) base++;
  }
  return base;
}

uint32_t lf_map(const FmIndex &fm, uint32_t row) {
  uint8_t c = fm.bwt[row];
  return fm.C[c] + occ(fm, c, row);
}

uint32_t row_to_pos(const FmIndex &fm, uint32_t row) {
  const uint32_t rate = fm.header.sa_sample_rate;
  uint32_t steps = 0;
  uint32_t cur = row;
  while (true) {
    if (cur % rate == 0) {
      uint32_t sample_idx = cur / rate;
      uint32_t sa_val = fm.sa_samples[sample_idx];
      return (sa_val + steps) % fm.text_len;
    }
    cur = lf_map(fm, cur);
    steps++;
  }
}

} // namespace

FmIndex build_fm_index(const std::vector<std::string> &protospacers,
                       const std::vector<FmIndex::Locus> &loci,
                       const FmIndexHeader &header) {
  if (protospacers.size() != loci.size()) {
    throw std::runtime_error("protospacers and loci size mismatch");
  }

  FmIndex fm;
  fm.header = header;

  // Build concatenated text with sentinel 4 after each protospacer to prevent cross-boundary matches.
  std::vector<uint8_t> text;
  text.reserve(protospacers.size() * (header.guide_length + 1));
  fm.loci.clear();
  fm.loci.reserve(protospacers.size());
  for (size_t i = 0; i < protospacers.size(); ++i) {
    for (char c : protospacers[i]) text.push_back(dna_to_sym(c));
    text.push_back(4); // sentinel
    fm.loci.push_back(loci[i]);
  }
  fm.text_len = static_cast<uint32_t>(text.size());

  // SA and BWT
  auto sa = build_sa(text);
  fm.bwt.resize(fm.text_len);
  for (uint32_t i = 0; i < fm.text_len; ++i) {
    uint32_t pos = sa[i];
    fm.bwt[i] = (pos == 0) ? 4 : text[pos - 1];
  }

  // C array
  std::array<uint32_t,5> counts{0,0,0,0,0};
  for (auto c : fm.bwt) counts[c]++;
  uint32_t cum = 0;
  for (size_t i = 0; i < counts.size(); ++i) {
    fm.C[i] = cum;
    cum += counts[i];
  }

  // Occ checkpoints
  uint32_t block = fm.header.occ_block;
  fm.occ_checkpoints.reserve((fm.text_len / block) + 2);
  std::array<uint32_t,5> running{0,0,0,0,0};
  for (uint32_t i = 0; i < fm.text_len; ++i) {
    if (i % block == 0) fm.occ_checkpoints.push_back(running);
    running[fm.bwt[i]]++;
  }
  fm.occ_checkpoints.push_back(running);

  // SA samples
  fm.sa_samples.reserve((fm.text_len / fm.header.sa_sample_rate) + 2);
  for (uint32_t i = 0; i < fm.text_len; ++i) {
    if (i % fm.header.sa_sample_rate == 0) fm.sa_samples.push_back(sa[i]);
  }

  return fm;
}

void save_fm_index(const FmIndex &fm, const std::string &path) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot open fm index for write: " + path);
  out.write(reinterpret_cast<const char *>(&fm.header), sizeof(FmIndexHeader));

  auto write_vec = [&out](const auto &v) {
    uint64_t n = v.size();
    out.write(reinterpret_cast<const char *>(&n), sizeof(n));
    out.write(reinterpret_cast<const char *>(v.data()), n * sizeof(typename std::decay_t<decltype(v)>::value_type));
  };

  write_vec(fm.bwt);
  write_vec(fm.occ_checkpoints);
  write_vec(fm.sa_samples);
  write_vec(fm.loci);
}

FmIndex load_fm_index(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open fm index for read: " + path);
  FmIndex fm;
  in.read(reinterpret_cast<char *>(&fm.header), sizeof(FmIndexHeader));

  auto read_vec = [&in](auto &v) {
    uint64_t n = 0; in.read(reinterpret_cast<char *>(&n), sizeof(n));
    using T = typename std::decay_t<decltype(v)>::value_type;
    v.resize(n);
    in.read(reinterpret_cast<char *>(v.data()), n * sizeof(T));
  };

  read_vec(fm.bwt);
  read_vec(fm.occ_checkpoints);
  read_vec(fm.sa_samples);
  read_vec(fm.loci);
  fm.text_len = static_cast<uint32_t>(fm.bwt.size());
  return fm;
}

std::vector<uint32_t> fm_search_exact_rows(const FmIndex &fm, const std::vector<uint8_t> &pattern) {
  uint32_t l = 0, r = fm.text_len;
  for (int i = static_cast<int>(pattern.size()) - 1; i >= 0; --i) {
    uint8_t c = pattern[i];
    l = fm.C[c] + occ(fm, c, l);
    r = fm.C[c] + occ(fm, c, r);
    if (l >= r) return {};
  }
  std::vector<uint32_t> rows;
  rows.reserve(r - l);
  for (uint32_t row = l; row < r; ++row) rows.push_back(row);
  return rows;
}

uint32_t fm_row_to_pos(const FmIndex &fm, uint32_t row) {
  return row_to_pos(fm, row);
}

namespace {

struct SearchState {
  int i;
  uint32_t l;
  uint32_t r;
  uint8_t mismatches;
};

} // namespace

std::vector<uint32_t> fm_search_hamming(const FmIndex &fm,
                                        const std::string &pattern,
                                        uint8_t max_mismatches) {
  std::vector<uint8_t> pat_enc(pattern.size());
  std::transform(pattern.begin(), pattern.end(), pat_enc.begin(), dna_to_sym);

  if (max_mismatches == 0) {
    return fm_search_exact_rows(fm, pat_enc);
  }

  std::vector<uint32_t> rows;
  std::vector<SearchState> stack;
  stack.reserve(1024);
  stack.push_back({static_cast<int>(pat_enc.size()) - 1, 0u, fm.text_len, 0u});

  while (!stack.empty()) {
    SearchState s = stack.back();
    stack.pop_back();

    if (s.i < 0) {
      for (uint32_t row = s.l; row < s.r; ++row) rows.push_back(row);
      continue;
    }

    for (uint8_t c = 0; c < 5; ++c) {
      uint8_t new_mm = s.mismatches + static_cast<uint8_t>(c != pat_enc[s.i]);
      if (new_mm > max_mismatches) continue;
      uint32_t l2 = fm.C[c] + occ(fm, c, s.l);
      uint32_t r2 = fm.C[c] + occ(fm, c, s.r);
      if (l2 >= r2) continue;
      stack.push_back({s.i - 1, l2, r2, new_mm});
    }
  }
  // dedup rows (could have duplicates when branching).
  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
  return rows;
}

} // namespace crispr_gpu
