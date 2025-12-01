#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/types.hpp"
#include "crispr_gpu/fm_index.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

namespace crispr_gpu {
namespace {

bool iupac_match(char base, char code) {
  base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
  code = static_cast<char>(std::toupper(static_cast<unsigned char>(code)));
  switch (code) {
    case 'N': return true;
    case 'A': return base == 'A';
    case 'C': return base == 'C';
    case 'G': return base == 'G';
    case 'T': return base == 'T';
    case 'R': return base == 'A' || base == 'G';
    case 'Y': return base == 'C' || base == 'T';
    case 'S': return base == 'G' || base == 'C';
    case 'W': return base == 'A' || base == 'T';
    case 'K': return base == 'G' || base == 'T';
    case 'M': return base == 'A' || base == 'C';
    case 'B': return base == 'C' || base == 'G' || base == 'T';
    case 'D': return base == 'A' || base == 'G' || base == 'T';
    case 'H': return base == 'A' || base == 'C' || base == 'T';
    case 'V': return base == 'A' || base == 'C' || base == 'G';
    default: return false;
  }
}

bool matches_pam(const std::string &segment, const std::string &pam) {
  if (segment.size() != pam.size()) return false;
  for (size_t i = 0; i < pam.size(); ++i) {
    if (!iupac_match(segment[i], pam[i])) return false;
  }
  return true;
}

struct FastaRecord {
  std::string name;
  std::string sequence;
};

std::vector<FastaRecord> read_fasta(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open FASTA: " + path);
  }
  std::vector<FastaRecord> records;
  std::string line;
  FastaRecord current;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') {
      if (!current.name.empty()) {
        records.push_back(current);
        current = FastaRecord{};
      }
      current.name = line.substr(1);
    } else {
      current.sequence += line;
    }
  }
  if (!current.name.empty()) records.push_back(current);
  return records;
}

void add_plus_strand_sites(const std::string &seq, uint32_t chrom_id, uint8_t guide_len,
                           const std::string &pam, std::vector<SiteRecord> &sites,
                           std::vector<std::string> *protospacers = nullptr,
                           std::vector<FmIndex::Locus> *loci = nullptr) {
  const size_t pam_len = pam.size();
  if (seq.size() < guide_len + pam_len) return;
  for (size_t i = 0; i + guide_len + pam_len <= seq.size(); ++i) {
    const std::string pam_seg = seq.substr(i + guide_len, pam_len);
    if (!matches_pam(pam_seg, pam)) continue;
    const std::string guide_seq = seq.substr(i, guide_len);
    try {
      uint64_t bits = encode_sequence_2bit(guide_seq);
      SiteRecord rec;
      rec.seq_bits = bits;
      rec.chrom_id = chrom_id;
      rec.pos = static_cast<uint32_t>(i);
      rec.strand = 0;
      sites.push_back(rec);
      if (protospacers && loci) {
        protospacers->push_back(guide_seq);
        loci->push_back({chrom_id, static_cast<uint32_t>(i), 0, static_cast<uint32_t>(sites.size() - 1)});
      }
    } catch (const std::runtime_error &) {
      continue; // skip sequences with ambiguous bases
    }
  }
}

void add_minus_strand_sites(const std::string &seq, uint32_t chrom_id, uint8_t guide_len,
                            const std::string &pam, std::vector<SiteRecord> &sites,
                            std::vector<std::string> *protospacers = nullptr,
                            std::vector<FmIndex::Locus> *loci = nullptr) {
  const std::string pam_rc = revcomp(pam);
  const size_t pam_len = pam.size();
  if (seq.size() < guide_len + pam_len) return;
  for (size_t i = 0; i + guide_len + pam_len <= seq.size(); ++i) {
    const std::string pam_seg = seq.substr(i, pam_len);
    if (!matches_pam(pam_seg, pam_rc)) continue;
    const std::string protospacer = seq.substr(i + pam_len, guide_len);
    const std::string protospacer_rc = revcomp(protospacer);
    try {
      uint64_t bits = encode_sequence_2bit(protospacer_rc);
      SiteRecord rec;
      rec.seq_bits = bits;
      rec.chrom_id = chrom_id;
      rec.pos = static_cast<uint32_t>(i + pam_len);
      rec.strand = 1;
      sites.push_back(rec);
      if (protospacers && loci) {
        protospacers->push_back(protospacer_rc);
        loci->push_back({chrom_id, static_cast<uint32_t>(i + pam_len), 1, static_cast<uint32_t>(sites.size() - 1)});
      }
    } catch (const std::runtime_error &) {
      continue;
    }
  }
}

} // namespace

GenomeIndex GenomeIndex::build(const std::string &fasta_path, const IndexParams &params) {
  if (params.guide_length == 0 || params.guide_length > 32) {
    throw std::runtime_error("guide_length must be between 1 and 32");
  }

  GenomeIndex out;
  out.meta_.guide_length = params.guide_length;
  out.meta_.pam = params.pam;
  out.meta_.both_strands = params.both_strands;

  const auto records = read_fasta(fasta_path);
  uint32_t chrom_id = 0;
  for (const auto &rec : records) {
    out.chroms_.push_back({rec.name, rec.sequence.size()});
    std::vector<std::string> protos;
    std::vector<FmIndex::Locus> loci;

    add_plus_strand_sites(rec.sequence, chrom_id, params.guide_length, params.pam, out.sites_, &protos, &loci);
    if (params.both_strands) {
      add_minus_strand_sites(rec.sequence, chrom_id, params.guide_length, params.pam, out.sites_, &protos, &loci);
    }

    FmIndexHeader hdr;
    hdr.guide_length = params.guide_length;
    hdr.pam = params.pam;
    hdr.both_strands = params.both_strands;
    hdr.occ_block = 128;
    hdr.sa_sample_rate = 32;
    auto fm = build_fm_index(protos, loci, hdr);
    out.fm_indices_.push_back(std::move(fm));
    ++chrom_id;
  }
  return out;
}

} // namespace crispr_gpu
