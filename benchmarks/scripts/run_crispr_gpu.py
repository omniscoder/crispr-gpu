#!/usr/bin/env python3
import argparse, json, subprocess, tempfile, time, pathlib, shutil

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--guides', required=True)
    ap.add_argument('--max-mm', type=int, default=4)
    ap.add_argument('--score-model', default='cfd')
    ap.add_argument('--backend', default='gpu')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    work = pathlib.Path(tempfile.mkdtemp(prefix='crispr_gpu_bench_'))
    idx = work / 'index.idx'

    t0 = time.time()
    subprocess.check_call(['crispr-gpu', 'index', '--fasta', args.fasta, '--pam', 'NGG', '--guide-length', '20', '--out', str(idx)])
    t1 = time.time()
    hits = work / 'hits.tsv'
    subprocess.check_call(['crispr-gpu', 'score', '--index', str(idx), '--guides', args.guides,
                           '--max-mm', str(args.max_mm), '--score-model', args.score_model,
                           '--backend', args.backend, '--output', str(hits)])
    t2 = time.time()

    result = {
        'tool': 'crispr-gpu',
        'backend': args.backend,
        'fasta': args.fasta,
        'guides': args.guides,
        'index_time_sec': t1 - t0,
        'score_time_sec': t2 - t1,
    }
    pathlib.Path(args.out).write_text(json.dumps(result, indent=2))
    shutil.rmtree(work)

if __name__ == '__main__':
    main()
