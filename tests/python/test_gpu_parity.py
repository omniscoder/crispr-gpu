import pytest
import crispr_gpu as cg


def build_toy_index(tmp_path):
    fasta = tmp_path / "toy.fa"
    fasta.write_text(">chr1\nAAAAGGGAAAA\n")
    params = cg.IndexParams()
    params.guide_length = 4
    params.pam = "NGG"
    params.both_strands = True
    return cg.GenomeIndex.build(str(fasta), params)


def sort_hits(hits):
    return sorted(
        hits,
        key=lambda h: (h.chrom_id, h.pos, h.strand, h.mismatches, round(h.score, 6)),
    )


@pytest.mark.skipif(not cg.cuda_available(), reason="CUDA backend not available")
def test_gpu_matches_cpu(tmp_path):
    idx = build_toy_index(tmp_path)
    guide = cg.Guide()
    guide.name = "g1"
    guide.sequence = "AAAA"
    guide.pam = "NGG"

    cpu_params = cg.EngineParams()
    cpu_params.max_mismatches = 4
    cpu_params.backend = cg.Backend.CPU

    gpu_params = cg.EngineParams()
    gpu_params.max_mismatches = 4
    gpu_params.backend = cg.Backend.GPU

    cpu_hits = sort_hits(cg.OffTargetEngine(idx, cpu_params).score_guide(guide))
    gpu_hits = sort_hits(cg.OffTargetEngine(idx, gpu_params).score_guide(guide))

    assert len(cpu_hits) == len(gpu_hits)
    for a, b in zip(cpu_hits, gpu_hits):
        assert (a.chrom_id, a.pos, a.strand, a.mismatches) == (
            b.chrom_id,
            b.pos,
            b.strand,
            b.mismatches,
        )
        assert pytest.approx(a.score, rel=1e-6, abs=1e-6) == b.score
