#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "crispr_gpu/genome_index.hpp"
#include "crispr_gpu/engine.hpp"
#include "crispr_gpu/version.hpp"

namespace py = pybind11;
using namespace crispr_gpu;

PYBIND11_MODULE(crispr_gpu, m) {
  m.attr("__version__") = CRISPR_GPU_VERSION;
  py::class_<Guide>(m, "Guide")
      .def(py::init<>())
      .def_readwrite("name", &Guide::name)
      .def_readwrite("sequence", &Guide::sequence)
      .def_readwrite("pam", &Guide::pam);

  py::enum_<ScoreModel>(m, "ScoreModel")
      .value("Hamming", ScoreModel::Hamming)
      .value("MIT", ScoreModel::MIT)
      .value("CFD", ScoreModel::CFD);

  py::enum_<Backend>(m, "Backend")
      .value("CPU", Backend::CPU)
      .value("GPU", Backend::GPU);

  py::class_<ScoreParams>(m, "ScoreParams")
      .def(py::init<>())
      .def_readwrite("model", &ScoreParams::model);

  py::class_<EngineParams>(m, "EngineParams")
      .def(py::init<>())
      .def_readwrite("max_mismatches", &EngineParams::max_mismatches)
      .def_readwrite("backend", &EngineParams::backend)
      .def_readwrite("score_params", &EngineParams::score_params);

  m.def("load_cfd_tables", &load_cfd_tables, "Load CFD/MM/position tables from JSON");

  py::class_<IndexParams>(m, "IndexParams")
      .def(py::init<>())
      .def_readwrite("guide_length", &IndexParams::guide_length)
      .def_readwrite("pam", &IndexParams::pam)
      .def_readwrite("both_strands", &IndexParams::both_strands);

  py::class_<OffTargetHit>(m, "OffTargetHit")
      .def_readonly("guide_name", &OffTargetHit::guide_name)
      .def_readonly("chrom_id", &OffTargetHit::chrom_id)
      .def_readonly("pos", &OffTargetHit::pos)
      .def_readonly("strand", &OffTargetHit::strand)
      .def_readonly("mismatches", &OffTargetHit::mismatches)
      .def_readonly("score", &OffTargetHit::score);

  py::class_<GenomeIndex>(m, "GenomeIndex")
      .def_static("load", &GenomeIndex::load)
      .def_static("build", &GenomeIndex::build)
      .def("save", &GenomeIndex::save)
      .def("num_sites", [](const GenomeIndex &g) { return g.sites().size(); })
      .def("guide_length", [](const GenomeIndex &g) { return g.meta().guide_length; })
      .def("pam", [](const GenomeIndex &g) { return g.meta().pam; });

  py::class_<OffTargetEngine>(m, "OffTargetEngine")
      .def(py::init<const GenomeIndex &, EngineParams>(), py::keep_alive<1, 2>())
      .def("score_guide", &OffTargetEngine::score_guide)
      .def("score_guides", &OffTargetEngine::score_guides);

  m.def("load_index", &GenomeIndex::load, "Load an index file");
  m.def("build_index", &GenomeIndex::build, "Build an index from FASTA",
        py::arg("fasta_path"), py::arg("params"));
  m.def("load_cfd_tables", &load_cfd_tables, "Load CFD/MM/position tables from JSON");
}
