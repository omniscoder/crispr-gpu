import json
import os
import subprocess
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _cli_binary() -> Path:
    # Matches existing convention in tests/python/test_api.py
    return _repo_root() / "build" / "crispr-gpu"


def test_cli_score_json_schema_and_sort(tmp_path):
    repo_root = _repo_root()
    binary = _cli_binary()
    assert binary.exists()

    fasta = tmp_path / "toy.fa"
    fasta.write_text(">chr1\nAAAAGGGAAAA\n")
    guides = tmp_path / "guides.tsv"
    # Include a header row; CLI should skip it.
    guides.write_text("name\tsequence\tpam\ng1\tAAAA\tNGG\ng2\tTTTT\tNGG\n")
    idx_path = tmp_path / "toy.idx"
    out_json = tmp_path / "hits.json"

    subprocess.check_call(
        [
            str(binary),
            "index",
            "--fasta",
            str(fasta),
            "--guide-length",
            "4",
            "--pam",
            "NGG",
            "--out",
            str(idx_path),
        ],
        cwd=repo_root,
    )

    subprocess.check_call(
        [
            str(binary),
            "score",
            "--index",
            str(idx_path),
            "--guides",
            str(guides),
            "--max-mm",
            "4",
            "--score-model",
            "hamming",
            "--backend",
            "cpu",
            "--search-backend",
            "brute",
            "--output-format",
            "json",
            "--sort",
            "--output",
            str(out_json),
        ],
        cwd=repo_root,
    )

    doc = json.loads(out_json.read_text())
    assert doc["schema"] == "crispr-gpu/score_result/v1"
    assert doc["schema_version"] == 1
    assert doc["tool"]["name"] == "crispr-gpu"
    assert doc["params"]["backend"] == "cpu"
    assert doc["params"]["sorted"] is True
    assert doc["index"]["guide_length"] == 4
    assert doc["input"]["num_guides"] == 2

    hits = doc["hits"]
    # Deterministic order: (guide, chrom, pos, strand, mismatches)
    keys = [(h["guide"], h["chrom"], h["pos"], h["strand"], h["mismatches"]) for h in hits]
    assert keys == sorted(keys)


def test_synthetic_bench_jsonl_schema(tmp_path):
    repo_root = _repo_root()
    script = repo_root / "benchmarks" / "run_synthetic.sh"
    jsonl = tmp_path / "synthetic.jsonl"

    env = dict(os.environ)
    env.update(
        {
            "BENCH_JSONL": str(jsonl),
            "BENCH_SCALE": "small",
            "GENOME_LEN": "10000",
            "GUIDE_COUNT": "1",
            "SKIP_GPU": "1",
            "CRISPR_GPU_WARMUP": "0",
        }
    )
    subprocess.check_call([str(script)], cwd=repo_root, env=env)

    lines = jsonl.read_text().strip().splitlines()
    assert len(lines) >= 1
    obj = json.loads(lines[0])
    assert obj["schema"] == "crispr-gpu/bench_run/v1"
    assert obj["schema_version"] == 1
    assert obj["tool"] == "crispr-gpu"
    assert obj["bench"] == "synthetic"
    assert obj["params"]["seed"] == 0
    assert obj["params"]["genome_len_bp"] == 10000
    assert obj["params"]["guides"] == 1

