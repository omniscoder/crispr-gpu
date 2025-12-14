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
#include <tuple>
#include <iomanip>

#include <nlohmann/json.hpp>

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
  std::string output_format{"tsv"}; // tsv|json|jsonl
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
  bool sort_hits{false};
};

static void print_usage() {
  std::cerr << "crispr-gpu " << CRISPR_GPU_VERSION << "\n";
  std::cerr << "Usage:\n";
  std::cerr << "  crispr-gpu index --fasta hg38.fa --pam NGG --guide-length 20 --out hg38.idx\n";
  std::cerr << "  crispr-gpu score --index hg38.idx --guides guides.tsv --max-mm 4 --score-model hamming --backend cpu|gpu --output hits.tsv\n";
  std::cerr << "  crispr-gpu score --search-backend brute|fmi  # fmi = exact K=0 only\n";
  std::cerr << "  crispr-gpu score --score-table table.json  # override MIT/CFD weights\n";
  std::cerr << "  crispr-gpu score --output-format tsv|json|jsonl [--sort]\n";
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

static std::string lower_ascii(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

static std::string score_model_to_string(ScoreModel m) {
  switch (m) {
    case ScoreModel::Hamming: return "hamming";
    case ScoreModel::MIT: return "mit";
    case ScoreModel::CFD: return "cfd";
  }
  return "unknown";
}

static std::string backend_to_string(Backend b) {
  switch (b) {
    case Backend::CPU: return "cpu";
    case Backend::GPU: return "gpu";
  }
  return "unknown";
}

static std::string search_backend_to_string(SearchBackend b) {
  switch (b) {
    case SearchBackend::BruteForce: return "brute";
    case SearchBackend::FMIndex: return "fmi";
  }
  return "unknown";
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
    // Skip common TSV header rows (e.g. "name sequence pam").
    {
      std::string n = lower_ascii(g.name);
      std::string s = lower_ascii(g.sequence);
      if ((n == "name" || n == "guide" || n == "id") &&
          (s == "sequence" || s == "seq")) {
        continue;
      }
    }
    if (g.pam.empty()) g.pam = "NGG";
    guides.push_back(g);
  }
  return guides;
}

static void sort_hits_in_place(std::vector<OffTargetHit> &hits) {
  std::sort(hits.begin(), hits.end(), [](const OffTargetHit &a, const OffTargetHit &b) {
    return std::tie(a.guide_name, a.chrom_id, a.pos, a.strand, a.mismatches) <
           std::tie(b.guide_name, b.chrom_id, b.pos, b.strand, b.mismatches);
  });
}

static std::string format_score_fixed(float score) {
  std::ostringstream ss;
  ss.setf(std::ios::fixed);
  ss << std::setprecision(8) << score;
  return ss.str();
}

static void write_hits_tsv(std::ostream &out,
                           const GenomeIndex &idx,
                           const std::vector<OffTargetHit> &hits) {
  out << "guide\tchrom\tpos\tstrand\tmismatches\tscore\n";
  const auto &chroms = idx.chromosomes();
  for (const auto &h : hits) {
    std::string chrom = (h.chrom_id < chroms.size()) ? chroms[h.chrom_id].name : std::to_string(h.chrom_id);
    out << h.guide_name << '\t' << chrom << '\t' << h.pos << '\t' << h.strand
        << '\t' << static_cast<int>(h.mismatches) << '\t' << format_score_fixed(h.score) << '\n';
  }
}

static void write_hits_jsonl(std::ostream &out,
                             const GenomeIndex &idx,
                             const std::vector<OffTargetHit> &hits) {
  // JSONL per-hit for streaming large results. Each line is a stable object (no outer wrapper).
  const auto &chroms = idx.chromosomes();
  for (const auto &h : hits) {
    nlohmann::json j;
    j["schema"] = "crispr-gpu/hit/v1";
    j["schema_version"] = 1;
    j["guide"] = h.guide_name;
    j["chrom"] = (h.chrom_id < chroms.size()) ? chroms[h.chrom_id].name : std::to_string(h.chrom_id);
    j["pos"] = h.pos;
    j["strand"] = std::string(1, h.strand);
    j["mismatches"] = static_cast<int>(h.mismatches);
    j["score"] = h.score;
    out << j.dump() << "\n";
  }
}

static void write_hits_json(std::ostream &out,
                            const CLIOptionsScore &opt,
                            const GenomeIndex &idx,
                            const std::vector<Guide> &guides,
                            const std::vector<OffTargetHit> &hits,
                            Backend effective_backend,
                            bool cuda_enabled) {
  nlohmann::json j;
  j["schema"] = "crispr-gpu/score_result/v1";
  j["schema_version"] = 1;
  j["tool"] = {
      {"name", "crispr-gpu"},
      {"version", std::string(CRISPR_GPU_VERSION)},
      {"cuda_enabled", cuda_enabled},
  };
  j["params"] = {
      {"backend", backend_to_string(effective_backend)},
      {"search_backend", search_backend_to_string(opt.search_backend)},
      {"score_model", score_model_to_string(opt.score_model)},
      {"max_mm", static_cast<int>(opt.max_mm)},
      {"score_table", opt.score_table.empty() ? nlohmann::json(nullptr) : nlohmann::json(opt.score_table)},
      {"sorted", opt.sort_hits},
  };
  j["index"] = {
      {"path", opt.index_path},
      {"guide_length", static_cast<int>(idx.meta().guide_length)},
      {"pam", idx.meta().pam},
      {"both_strands", idx.meta().both_strands},
      {"num_sites", static_cast<std::uint64_t>(idx.sites().size())},
      {"num_chromosomes", static_cast<std::uint64_t>(idx.chromosomes().size())},
  };
  j["input"] = {
      {"guides_path", opt.guides_path},
      {"num_guides", static_cast<std::uint64_t>(guides.size())},
  };

  const auto &chroms = idx.chromosomes();
  nlohmann::json arr = nlohmann::json::array();
  for (const auto &h : hits) {
    nlohmann::json row;
    row["guide"] = h.guide_name;
    row["chrom"] = (h.chrom_id < chroms.size()) ? chroms[h.chrom_id].name : std::to_string(h.chrom_id);
    row["pos"] = h.pos;
    row["strand"] = std::string(1, h.strand);
    row["mismatches"] = static_cast<int>(h.mismatches);
    row["score"] = h.score;
    arr.push_back(std::move(row));
  }
  j["hits"] = std::move(arr);
  j["summary"] = {
      {"num_hits", static_cast<std::uint64_t>(hits.size())},
  };

  out << j.dump(2, ' ', true) << "\n";
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
        if (arg == "--output-format" && i + 1 < argc) { opt.output_format = lower_ascii(argv[++i]); continue; }
        if (arg == "--score-table" && i + 1 < argc) { opt.score_table = argv[++i]; continue; }
        if (arg == "--max-mm" && i + 1 < argc) { opt.max_mm = static_cast<uint8_t>(std::stoi(argv[++i])); continue; }
        if (arg == "--score-model" && i + 1 < argc) { opt.score_model = parse_score_model(argv[++i]); continue; }
        if (arg == "--backend" && i + 1 < argc) { opt.backend = parse_backend(argv[++i]); continue; }
        if (arg == "--search-backend" && i + 1 < argc) { opt.search_backend = parse_search_backend(argv[++i]); continue; }
        if (arg == "--sort") { opt.sort_hits = true; continue; }
      }
      if (opt.index_path.empty() || opt.guides_path.empty()) {
        throw std::runtime_error("--index and --guides are required");
      }
      if (opt.output_format != "tsv" && opt.output_format != "json" && opt.output_format != "jsonl") {
        throw std::runtime_error("--output-format must be one of: tsv, json, jsonl");
      }
      const bool cuda_enabled =
#ifdef CRISPR_GPU_ENABLE_CUDA
          true
#else
          false
#endif
          ;

      Backend effective_backend = opt.backend;
      if (effective_backend == Backend::GPU && !cuda_available()) {
        std::cerr << "CUDA backend not available at runtime; falling back to CPU.\n";
        effective_backend = Backend::CPU;
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
      ep.backend = effective_backend;
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
      if (opt.sort_hits) {
        ScopedTimer t_sort("cli.score.sort_hits", timing_enabled());
        sort_hits_in_place(hits);
      }

      std::ostream *out = &std::cout;
      std::ofstream outfile;
      if (!opt.output_path.empty()) {
        outfile.open(opt.output_path);
        out = &outfile;
      }
      {
        ScopedTimer t_write("cli.score.write_output", timing_enabled());
        if (opt.output_format == "tsv") {
          write_hits_tsv(*out, idx, hits);
        } else if (opt.output_format == "jsonl") {
          write_hits_jsonl(*out, idx, hits);
        } else {
          write_hits_json(*out, opt, idx, guides, hits, effective_backend, cuda_enabled);
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
