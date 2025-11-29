import os
import tempfile
import subprocess
import pytest
from pathlib import Path
import crispr_gpu as cg


def build_toy_index(tmp_path):
    fasta = tmp_path / "toy.fa"
    fasta.write_text(">chr1\nAAAAGGGAAAA\n")
    params = cg.IndexParams()
    params.guide_length = 4
    params.pam = "NGG"
    params.both_strands = True
    idx = cg.GenomeIndex.build(str(fasta), params)
    return idx


def test_python_import():
    import crispr_gpu  # noqa: F401


def test_python_build_and_score_cpu(tmp_path):
    idx = build_toy_index(tmp_path)
    g = cg.Guide()
    g.name = "g1"
    g.sequence = "AAAA"
    g.pam = "NGG"

    params = cg.EngineParams()
    params.max_mismatches = 4
    params.backend = cg.Backend.CPU

    eng = cg.OffTargetEngine(idx, params)
    hits = eng.score_guide(g)
    assert len(hits) >= 1
    assert hits[0].mismatches == 0


def test_cli_cpu_end_to_end(tmp_path):
    fasta = tmp_path / "toy.fa"
    fasta.write_text(">chr1\nAAAAGGGAAAA\n")
    guides = tmp_path / "guides.tsv"
    guides.write_text("g1\tAAAA\tNGG\n")
    idx_path = tmp_path / "toy.idx"
    hits_path = tmp_path / "hits.tsv"

    # Assume build dir relative to repo root (tests/python/../..)
    repo_root = Path(__file__).resolve().parents[2]
    binary = repo_root / "build" / "crispr-gpu"
    subprocess.check_call([
        str(binary),
        "index", "--fasta", str(fasta), "--guide-length", "4", "--pam", "NGG", "--out", str(idx_path)
    ])
    subprocess.check_call([
        str(binary),
        "score", "--index", str(idx_path), "--guides", str(guides),
        "--max-mm", "4", "--score-model", "hamming", "--backend", "cpu", "--output", str(hits_path)
    ])
    assert hits_path.exists()
    content = hits_path.read_text().strip().splitlines()
    assert len(content) >= 2  # header + at least one hit
