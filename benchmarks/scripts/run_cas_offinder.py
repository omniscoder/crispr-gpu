#!/usr/bin/env python3
import argparse, json, subprocess, tempfile, pathlib, time, shutil

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cas-offinder', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--guides', required=True)
    ap.add_argument('--pam', default='NGG')
    ap.add_argument('--max-mm', type=int, default=4)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    work = pathlib.Path(tempfile.mkdtemp(prefix='cas_offinder_bench_'))
    # Cas-OFFinder input format: pattern, genome, output
    input_txt = work / 'input.txt'
    # pattern: N{guide_len}[pam], here 20 + pam
    pattern = 'N' * 20 + args.pam
    input_txt.write_text('\n'.join([pattern, args.fasta, args.guides]))

    t0 = time.time()
    subprocess.check_call([args.cas_offinder, str(input_txt), 'C', str(args.max_mm), str(work/'hits.txt')])
    t1 = time.time()

    result = {
        'tool': 'cas-offinder',
        'backend': 'gpu',
        'fasta': args.fasta,
        'guides': args.guides,
        'score_time_sec': t1 - t0,
    }
    pathlib.Path(args.out).write_text(json.dumps(result, indent=2))
    shutil.rmtree(work)

if __name__ == '__main__':
    main()
