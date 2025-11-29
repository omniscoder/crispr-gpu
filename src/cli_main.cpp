#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/engine.hpp"
#include "crispr_gpu/version.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>

using namespace crispr_gpu;

struct CLIOptionsIndex {
  std::string fasta;
  std::string pam{"NGG"};
  uint8_t guide_len{20};
  bool both_strands{true};
  std::string out{"index.idx"};
};

struct CLIOptionsScore {
  std::string index_path;
  std::string guides_path;
  std::string output_path{""};
  std::string cfd_table{""};
  uint8_t max_mm{4};
  ScoreModel score_model{ScoreModel::Hamming};
  Backend backend{
#ifdef CRISPR_GPU_ENABLE_CUDA
      Backend::GPU
#else
      Backend::CPU
#endif
  };
};

static void print_usage() {
  std::cerr << "crispr-gpu " << CRISPR_GPU_VERSION << "\n";
  std::cerr << "Usage:\n";
  std::cerr << "  crispr-gpu index --fasta hg38.fa --pam NGG --guide-length 20 --out hg38.idx\n";
  std::cerr << "  crispr-gpu score --index hg38.idx --guides guides.tsv --max-mm 4 --score-model hamming --backend cpu|gpu --output hits.tsv\n";
  std::cerr << "  crispr-gpu score --cfd-table cfd.json  # override CFD weights\n";
  std::cerr << "  crispr-gpu --version\n";
}

static ScoreModel parse_score_model(const std::string &s) {
  std::string l = s;
  std::transform(l.begin(), l.end(), l.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (l == "hamming") return ScoreModel::Hamming;
  if (l == "mit") return ScoreModel::MIT;
  if (l == "cfd") return ScoreModel::CFD;
  throw std::runtime_error("Unknown score model: " + s);
}

static Backend parse_backend(const std::string &s) {
  std::string l = s;
  std::transform(l.begin(), l.end(), l.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (l == "cpu") return Backend::CPU;
  if (l == "gpu") return Backend::GPU;
  throw std::runtime_error("Unknown backend: " + s);
}

static std::vector<Guide> read_guides(const std::string &path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Failed to open guides file: " + path);
  std::vector<Guide> guides;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    Guide g;
    ss >> g.name >> g.sequence >> g.pam;
    if (g.name.empty() || g.sequence.empty()) continue;
    if (g.pam.empty()) g.pam = "NGG";
    guides.push_back(g);
  }
  return guides;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    print_usage();
    return 1;
  }

  std::string sub = argv[1];
  try {
    if (sub == "--version" || sub == "-V") {
      std::cout << CRISPR_GPU_VERSION << "\n";
      return 0;
    }
    if (sub == "index") {
      CLIOptionsIndex opt;
      for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--fasta" && i + 1 < argc) { opt.fasta = argv[++i]; continue; }
        if (arg == "--pam" && i + 1 < argc) { opt.pam = argv[++i]; continue; }
        if (arg == "--guide-length" && i + 1 < argc) { opt.guide_len = static_cast<uint8_t>(std::stoi(argv[++i])); continue; }
        if (arg == "--out" && i + 1 < argc) { opt.out = argv[++i]; continue; }
        if (arg == "--both-strands") { opt.both_strands = true; continue; }
        if (arg == "--plus-only") { opt.both_strands = false; continue; }
      }
      if (opt.fasta.empty()) {
        throw std::runtime_error("--fasta is required");
      }
      IndexParams params;
      params.guide_length = opt.guide_len;
      params.pam = opt.pam;
      params.both_strands = opt.both_strands;
      auto idx = GenomeIndex::build(opt.fasta, params);
      idx.save(opt.out);
      std::cerr << "Index written to " << opt.out << " with " << idx.sites().size() << " sites\n";
      return 0;
    } else if (sub == "score") {
      CLIOptionsScore opt;
      for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--index" && i + 1 < argc) { opt.index_path = argv[++i]; continue; }
        if (arg == "--guides" && i + 1 < argc) { opt.guides_path = argv[++i]; continue; }
        if (arg == "--output" && i + 1 < argc) { opt.output_path = argv[++i]; continue; }
        if (arg == "--cfd-table" && i + 1 < argc) { opt.cfd_table = argv[++i]; continue; }
        if (arg == "--max-mm" && i + 1 < argc) { opt.max_mm = static_cast<uint8_t>(std::stoi(argv[++i])); continue; }
        if (arg == "--score-model" && i + 1 < argc) { opt.score_model = parse_score_model(argv[++i]); continue; }
        if (arg == "--backend" && i + 1 < argc) { opt.backend = parse_backend(argv[++i]); continue; }
      }
      if (opt.index_path.empty() || opt.guides_path.empty()) {
        throw std::runtime_error("--index and --guides are required");
      }
      auto idx = GenomeIndex::load(opt.index_path);
      auto guides = read_guides(opt.guides_path);
      EngineParams ep;
      ep.max_mismatches = opt.max_mm;
      ep.score_params.model = opt.score_model;
      ep.backend = opt.backend;
      if (!opt.cfd_table.empty()) {
        load_cfd_tables(opt.cfd_table);
      }
      OffTargetEngine engine(idx, ep);
      auto hits = engine.score_guides(guides);

      std::ostream *out = &std::cout;
      std::ofstream outfile;
      if (!opt.output_path.empty()) {
        outfile.open(opt.output_path);
        out = &outfile;
      }
      *out << "guide\tchrom\tpos\tstrand\tmismatches\tscore\n";
      for (const auto &h : hits) {
        const auto &chroms = idx.chromosomes();
        std::string chrom = (h.chrom_id < chroms.size()) ? chroms[h.chrom_id].name : std::to_string(h.chrom_id);
        *out << h.guide_name << '\t' << chrom << '\t' << h.pos << '\t' << h.strand
             << '\t' << static_cast<int>(h.mismatches) << '\t' << h.score << '\n';
      }
      return 0;
    } else {
      print_usage();
      return 1;
    }
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
