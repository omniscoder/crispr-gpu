#!/usr/bin/env python3
import argparse, json, glob
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', nargs='*', default=None, help='JSON result files')
    args = ap.parse_args()

    files = args.inputs or glob.glob('*.json')
    rows = []
    for f in files:
        with open(f) as fh:
            rows.append(json.load(fh))
    df = pd.DataFrame(rows)
    print(df[['tool','backend','guides','score_time_sec','index_time_sec']])
    if not df.empty:
        print('\nRuntime vs guides:')
        print(df.groupby(['tool'])['score_time_sec'].describe())

if __name__ == '__main__':
    main()
