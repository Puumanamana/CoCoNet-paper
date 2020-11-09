#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import h5py


def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata', type=str)
    parser.add_argument('--h5', type=str)
    parser.add_argument('--csv', action='store_true')
    args = parser.parse_args()

    return args

def summarize_metadata(path, csv=False):
    meta = pd.read_csv(path)

    name = Path(path).parent.name
    n_viruses = len(meta.V_id.unique())
    n_ctg = len(meta)
    n_ctg_gt1k = sum(meta['size'] >= 1024)    
    n_ctg_gt2k = sum(meta['size'] >= 2048)
    n_frag_virus = sum(meta.C_id.str.endswith('|1'))

    if csv:
        params = re.split(r'[\-_]', name)
        params = [params[i] for i in [1, 3, 5, 6]]
        entry = params + [n_viruses, n_ctg, n_ctg_gt1k, n_ctg_gt2k, n_frag_virus]
        print(','.join(map(str, entry)))
        return

    info = [
        ('=== simulation ===', name),
        ('n_viruses', f'{n_viruses:,} (fragmented: {n_frag_virus/n_viruses:.1%})'),
        ('n_contigs', f'{n_ctg:,} (>1024: {n_ctg_gt1k/n_ctg:.1%} >2048: {n_ctg_gt2k/n_ctg:.1%})'),
        ('fragmentation', f'{n_ctg/n_viruses:.2f}'),
    ]

    print('\n'.join(f'{n:>25}: {v}' for (n, v) in info))

def summarize_abundance(path, csv=False):
    with h5py.File(path, 'r') as handle:
        n_contigs = len(handle)
        n_samples = next(iter(handle.values())).shape[0]

        coverage = np.zeros((n_contigs, n_samples))
        ctg_sizes = np.zeros(n_contigs)

        for i, cov in enumerate(handle.values()):
            coverage[i] = cov[:].sum(axis=1)
            ctg_sizes[i] = cov.shape[1]

    xcoverage = coverage.sum(axis=0) / ctg_sizes[:, None].sum()
    gt2kb = ctg_sizes >= 2048
    xcov_per_ctg_gt2kb = (coverage[gt2kb, :] / ctg_sizes[gt2kb, None]).mean(axis=0)

    sim = re.split(r'[_\-]', Path(path).parent.name)

    if csv:
        entry = [sim[1], sim[3], sim[5], sim[6], n_contigs, sum(gt2kb),
                 xcoverage.mean(), xcov_per_ctg_gt2kb.mean()]
        print(','.join(map(str, entry)))
    else:
        print(f'======= Simulation: {n_samples} samples ({sim[5]}) =======')
        print(f'Mean overall coverage: {xcoverage.mean().round(4)}')
        print(f'Mean contig coverage (> 2kb): {xcov_per_ctg_gt2kb.mean().round(4)}')
        print(f'#contigs: (>0) {n_contigs:,} (>2048) {sum(gt2kb):,}')

def main():

    args = parse_args()

    if args.metadata is not None:
        summarize_metadata(args.metadata, csv=args.csv)
    if args.h5 is not None:
        summarize_abundance(args.h5, csv=args.csv)
    
        
if __name__ == '__main__':
    main()
