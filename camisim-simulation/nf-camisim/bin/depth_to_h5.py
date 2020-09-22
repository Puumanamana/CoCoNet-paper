import argparse
from pathlib import Path

import pandas as pd
import h5py


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata', type=str)
    parser.add_argument('--n-samples', type=int)    
    args = parser.parse_args()

    return args

def main():

    args = parse_args()

    metadata = pd.read_csv(args.metadata)

    info = metadata.groupby("V_id").agg(list)
    sizes = metadata.groupby("V_id")['size'].first()

    cov_vir_h5 = h5py.File("coverage_virus.h5","w")
    cov_ctg_h5 = h5py.File("coverage_contigs.h5","w")

    for filename in Path('.').glob('*.txt'):
        virus = filename.stem

        if virus not in info.index:
            continue

        coverage = pd.read_csv(
            filename, sep='\\t', header=None,
            usecols=range(1, args.n_samples+2), dtype=int
        ).set_index(1)
        coverage.index -= 1 # since samtools positions are 1-based
        coverage.columns = [f'sample_{i+1}' for i in range(args.n_samples)]
        coverage = (coverage
            .reindex(index=range(sizes.loc[virus]))
            .fillna(0)
            .to_numpy())

        cov_vir_h5.create_dataset(virus, data=coverage)

        for ctg, start, end in zip(*info.loc[virus,["C_id","start","end"]]):
            cov_ctg_h5.create_dataset(ctg, data=coverage[:,start:end+1])

    cov_ctg_h5.close()
    cov_vir_h5.close()

if __name__ == '__main__':
    main()
