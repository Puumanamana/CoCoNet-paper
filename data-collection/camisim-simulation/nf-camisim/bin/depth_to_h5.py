#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd
import h5py


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata', type=str)
    parser.add_argument('--genome-sizes', type=str)    
    parser.add_argument('--n-samples', type=int)    
    args = parser.parse_args()

    return args

def main():

    args = parse_args()

    metadata = pd.read_csv(args.metadata)

    info = metadata.groupby("V_id").agg(list)
    sizes = pd.read_csv(args.genome_sizes, index_col=0).iloc[:, 0]

    cov_vir_h5 = h5py.File("coverage_virus.h5","w")
    cov_ctg_h5 = h5py.File("coverage_contigs.h5","w")

    for filename in Path('.').glob('*.txt'):
        virus = filename.stem

        if virus not in info.index:
            continue

        coverage = pd.read_csv(
            filename, sep='\t', header=None, dtype=int,
            usecols=range(1, args.n_samples+2),
            names=['pos']+[f'sample-{i+1}' for i in range(args.n_samples)],
        ).set_index('pos')

        coverage = (coverage
                    .reindex(index=range(1, 1+sizes.loc[virus]))
                    .fillna(0)
                    .to_numpy().T)

        cov_vir_h5.create_dataset(virus, data=coverage)

        for ctg, start, end in zip(*info.loc[virus, ["C_id", "start", "end"]]):
            cov_ctg_h5.create_dataset(ctg, data=coverage[:,start-1:end])

    cov_ctg_h5.close()
    cov_vir_h5.close()

if __name__ == '__main__':
    main()
