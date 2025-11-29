#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/types.hpp"

#include <fstream>
#include <stdexcept>
#include <cstring>

namespace crispr_gpu {
namespace {
struct IndexHeader {
  char magic[8];
  uint32_t version;
  uint8_t guide_length;
  uint8_t pam_length;
  uint8_t both_strands;
  uint8_t reserved;
  uint64_t num_chroms;
  uint64_t num_sites;
};

constexpr const char *kMagic = "CRSPRDX"; // 7 chars + null padding
constexpr uint32_t kVersion = 1;
}

GenomeIndex GenomeIndex::load(const std::string &index_path) {
  std::ifstream in(index_path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("Failed to open index: " + index_path);
  }

  IndexHeader hdr{};
  in.read(reinterpret_cast<char *>(&hdr), sizeof(hdr));
  if (!in) throw std::runtime_error("Index truncated: header");

  if (std::strncmp(hdr.magic, kMagic, 7) != 0) {
    throw std::runtime_error("Invalid index magic");
  }
  if (hdr.version != kVersion) {
    throw std::runtime_error("Unsupported index version");
  }

  GenomeIndex idx;
  idx.meta_.guide_length = hdr.guide_length;
  idx.meta_.both_strands = hdr.both_strands != 0;

  std::string pam(hdr.pam_length, '\0');
  in.read(&pam[0], hdr.pam_length);
  if (!in) throw std::runtime_error("Index truncated: pam");
  idx.meta_.pam = pam;

  idx.chroms_.reserve(hdr.num_chroms);
  for (uint64_t i = 0; i < hdr.num_chroms; ++i) {
    uint16_t name_len = 0;
    in.read(reinterpret_cast<char *>(&name_len), sizeof(name_len));
    if (!in) throw std::runtime_error("Index truncated: chrom name length");
    std::string name(name_len, '\0');
    in.read(&name[0], name_len);
    if (!in) throw std::runtime_error("Index truncated: chrom name");
    uint64_t chrom_len = 0;
    in.read(reinterpret_cast<char *>(&chrom_len), sizeof(chrom_len));
    if (!in) throw std::runtime_error("Index truncated: chrom length");
    idx.chroms_.push_back({name, chrom_len});
  }

  idx.sites_.reserve(hdr.num_sites);
  for (uint64_t i = 0; i < hdr.num_sites; ++i) {
    SiteRecord rec{};
    in.read(reinterpret_cast<char *>(&rec.seq_bits), sizeof(rec.seq_bits));
    in.read(reinterpret_cast<char *>(&rec.chrom_id), sizeof(rec.chrom_id));
    in.read(reinterpret_cast<char *>(&rec.pos), sizeof(rec.pos));
    in.read(reinterpret_cast<char *>(&rec.strand), sizeof(rec.strand));
    if (!in) throw std::runtime_error("Index truncated: sites");
    idx.sites_.push_back(rec);
  }

  return idx;
}

void GenomeIndex::save(const std::string &index_path) const {
  std::ofstream out(index_path, std::ios::binary);
  if (!out) {
    throw std::runtime_error("Failed to write index: " + index_path);
  }

  IndexHeader hdr{};
  std::memset(&hdr, 0, sizeof(hdr));
  std::memcpy(hdr.magic, kMagic, 7);
  hdr.version = kVersion;
  hdr.guide_length = meta_.guide_length;
  hdr.pam_length = static_cast<uint8_t>(meta_.pam.size());
  hdr.both_strands = meta_.both_strands ? 1 : 0;
  hdr.num_chroms = static_cast<uint64_t>(chroms_.size());
  hdr.num_sites = static_cast<uint64_t>(sites_.size());

  out.write(reinterpret_cast<const char *>(&hdr), sizeof(hdr));
  out.write(meta_.pam.data(), meta_.pam.size());

  for (const auto &chrom : chroms_) {
    uint16_t name_len = static_cast<uint16_t>(chrom.name.size());
    out.write(reinterpret_cast<const char *>(&name_len), sizeof(name_len));
    out.write(chrom.name.data(), name_len);
    out.write(reinterpret_cast<const char *>(&chrom.length), sizeof(chrom.length));
  }

  for (const auto &rec : sites_) {
    out.write(reinterpret_cast<const char *>(&rec.seq_bits), sizeof(rec.seq_bits));
    out.write(reinterpret_cast<const char *>(&rec.chrom_id), sizeof(rec.chrom_id));
    out.write(reinterpret_cast<const char *>(&rec.pos), sizeof(rec.pos));
    out.write(reinterpret_cast<const char *>(&rec.strand), sizeof(rec.strand));
  }

  if (!out) {
    throw std::runtime_error("Failed while writing index: " + index_path);
  }
}

} // namespace crispr_gpu
