#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('metadata', type=str)
    args = parser.parse_args()

    return args

def main():
    '''
    '''

    args = parse_args()

    meta = pd.read_csv(args.metadata)

    name = Path(args.metadata).parent.name
    n_viruses = len(meta.V_id.unique())
    n_ctg = len(meta)
    n_ctg_gt2k = sum(meta['size'] >= 2048)
    n_frag_virus = sum(meta.C_id.str.endswith('|1'))

    info = [
        ('=== simulation ===', name),
        ('n_viruses', f'{n_viruses:,} (fragmented: {n_frag_virus/n_viruses:.1%})'),
        ('n_contigs', f'{n_ctg:,} (>2048: {n_ctg_gt2k/n_ctg:.1%})'),
        ('fragmentation', f'{n_ctg/n_viruses:.2f}'),
    ]

    print('\n'.join(f'{n:>25}: {v}' for (n, v) in info))

if __name__ == '__main__':
    main()
