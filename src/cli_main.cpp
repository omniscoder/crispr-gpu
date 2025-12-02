#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/engine.hpp"
#include "crispr_gpu/version.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>

using namespace crispr_gpu;

static bool timing_enabled() {
  static int cached = -1;
  if (cached == -1) {
    const char *env = std::getenv("CRISPR_GPU_TIMING");
    cached = (env && env[0] != '0') ? 1 : 0;
  }
  return cached == 1;
}

struct ScopedTimer {
  const char *name;
  bool active;
  std::chrono::steady_clock::time_point start;
  ScopedTimer(const char *n, bool on) : name(n), active(on), start(std::chrono::steady_clock::now()) {}
  ~ScopedTimer() {
    if (!active) return;
    auto end = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::fprintf(stderr, "[timing] %s: %.3f ms\n", name, ms);
  }
};

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
  std::string score_table{""};
  uint8_t max_mm{4};
  ScoreModel score_model{ScoreModel::Hamming};
  Backend backend{
#ifdef CRISPR_GPU_ENABLE_CUDA
      Backend::GPU
#else
      Backend::CPU
#endif
  };
  SearchBackend search_backend{SearchBackend::BruteForce};
};

static void print_usage() {
  std::cerr << "crispr-gpu " << CRISPR_GPU_VERSION << "\n";
  std::cerr << "Usage:\n";
  std::cerr << "  crispr-gpu index --fasta hg38.fa --pam NGG --guide-length 20 --out hg38.idx\n";
  std::cerr << "  crispr-gpu score --index hg38.idx --guides guides.tsv --max-mm 4 --score-model hamming --backend cpu|gpu --output hits.tsv\n";
  std::cerr << "  crispr-gpu score --search-backend brute|fmi  # fmi = exact K=0 only\n";
  std::cerr << "  crispr-gpu score --score-table table.json  # override MIT/CFD weights\n";
  std::cerr << "  crispr-gpu warmup  # warm CUDA context (no-op if CUDA disabled)\n";
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

static SearchBackend parse_search_backend(const std::string &s) {
  std::string l = s;
  std::transform(l.begin(), l.end(), l.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (l == "brute" || l == "bruteforce" || l == "scan") return SearchBackend::BruteForce;
  if (l == "fmi" || l == "fmindex") return SearchBackend::FMIndex;
  throw std::runtime_error("Unknown search backend: " + s);
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
    if (sub == "warmup") {
#ifdef CRISPR_GPU_ENABLE_CUDA
      cuda_warmup();
      std::cerr << "CUDA warmup done.\n";
#else
      std::cerr << "CUDA not enabled; warmup is a no-op.\n";
#endif
      return 0;
    }
    if (sub == "index") {
      ScopedTimer t_total("cli.index.total", timing_enabled());
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
      GenomeIndex idx;
      {
        ScopedTimer t_build("cli.index.build", timing_enabled());
        idx = GenomeIndex::build(opt.fasta, params);
      }
      {
        ScopedTimer t_save("cli.index.save", timing_enabled());
        idx.save(opt.out);
      }
      std::cerr << "Index written to " << opt.out << " with " << idx.sites().size() << " sites\n";
      return 0;
    } else if (sub == "score") {
      ScopedTimer t_total("cli.score.total", timing_enabled());
      CLIOptionsScore opt;
      for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--index" && i + 1 < argc) { opt.index_path = argv[++i]; continue; }
        if (arg == "--guides" && i + 1 < argc) { opt.guides_path = argv[++i]; continue; }
        if (arg == "--output" && i + 1 < argc) { opt.output_path = argv[++i]; continue; }
        if (arg == "--score-table" && i + 1 < argc) { opt.score_table = argv[++i]; continue; }
        if (arg == "--max-mm" && i + 1 < argc) { opt.max_mm = static_cast<uint8_t>(std::stoi(argv[++i])); continue; }
        if (arg == "--score-model" && i + 1 < argc) { opt.score_model = parse_score_model(argv[++i]); continue; }
        if (arg == "--backend" && i + 1 < argc) { opt.backend = parse_backend(argv[++i]); continue; }
        if (arg == "--search-backend" && i + 1 < argc) { opt.search_backend = parse_search_backend(argv[++i]); continue; }
      }
      if (opt.index_path.empty() || opt.guides_path.empty()) {
        throw std::runtime_error("--index and --guides are required");
      }
      GenomeIndex idx;
      {
        ScopedTimer t_load_index("cli.score.load_index", timing_enabled());
        idx = GenomeIndex::load(opt.index_path);
      }
      std::vector<Guide> guides;
      {
        ScopedTimer t_read_guides("cli.score.read_guides", timing_enabled());
        guides = read_guides(opt.guides_path);
      }
      EngineParams ep;
      ep.max_mismatches = opt.max_mm;
      ep.score_params.model = opt.score_model;
      ep.score_params.table_path = opt.score_table;
      ep.backend = opt.backend;
      ep.search_backend = opt.search_backend;
      if (!opt.score_table.empty()) {
        ScopedTimer t_cfd("cli.score.load_table", timing_enabled());
        load_cfd_tables(opt.score_table);
      }
      OffTargetEngine engine(idx, ep);
      std::vector<OffTargetHit> hits;
      {
        ScopedTimer t_score("cli.score.score_guides", timing_enabled());
        hits = engine.score_guides(guides);
      }

      std::ostream *out = &std::cout;
      std::ofstream outfile;
      if (!opt.output_path.empty()) {
        outfile.open(opt.output_path);
        out = &outfile;
      }
      {
        ScopedTimer t_write("cli.score.write_output", timing_enabled());
        *out << "guide\tchrom\tpos\tstrand\tmismatches\tscore\n";
        for (const auto &h : hits) {
          const auto &chroms = idx.chromosomes();
          std::string chrom = (h.chrom_id < chroms.size()) ? chroms[h.chrom_id].name : std::to_string(h.chrom_id);
          *out << h.guide_name << '\t' << chrom << '\t' << h.pos << '\t' << h.strand
               << '\t' << static_cast<int>(h.mismatches) << '\t' << h.score << '\n';
        }
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
