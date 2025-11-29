#!/usr/bin/env python3
import argparse, json, subprocess, tempfile, pathlib, time, shutil

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--crispritz', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--guides', required=True)
    ap.add_argument('--pam', default='NGG')
    ap.add_argument('--max-mm', type=int, default=4)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    work = pathlib.Path(tempfile.mkdtemp(prefix='crispritz_bench_'))
    # Build index
    idx_dir = work / 'idx'
    t0 = time.time()
    subprocess.check_call([args.crispritz, 'index', '-g', args.fasta, '-o', str(idx_dir), '-p', args.pam])
    t1 = time.time()
    # Score
    t2 = t1
    try:
      subprocess.check_call([args.crispritz, 'search', '-i', str(idx_dir), '-g', args.guides, '-m', str(args.max_mm), '-o', str(work/'hits')])
      t2 = time.time()
    finally:
      result = {
          'tool': 'crispritz',
          'backend': 'gpu',
          'fasta': args.fasta,
          'guides': args.guides,
          'index_time_sec': t1 - t0,
          'score_time_sec': t2 - t1 if t2>t1 else None,
      }
      pathlib.Path(args.out).write_text(json.dumps(result, indent=2))
      shutil.rmtree(work)

if __name__ == '__main__':
    main()
