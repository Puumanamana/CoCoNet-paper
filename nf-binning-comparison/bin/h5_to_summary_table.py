#!/usr/bin/env python3

from pathlib import Path
import argparse
import h5py
import numpy as np
import pandas as pd


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--abundance', type=str)
    parser.add_argument('--fasta', type=str)        
    parser.add_argument('--suffix', type=str)    
    args = parser.parse_args()

    return args

def h5_to_metabat2(h5_path):

    h5data = h5py.File(h5_path, 'r')

    ctg_len  = pd.Series({ctg: h5data[ctg].shape[1] for ctg in h5data.keys()}, name='ctg_len')
    contigs = ctg_len.index
    n_samples = h5data[contigs[0]].shape[0]

    coverage = np.zeros([len(contigs), 1+2*n_samples])

    for i, ctg in enumerate(contigs):
        data = h5data.get(ctg)[:]
        mu, sig2 = (np.mean(data,axis=1), np.var(data,axis=1))
        coverage[i, 0] = np.sum(data) / data.shape[1]
        coverage[i, 1::2] = mu
        coverage[i, 2::2] = sig2

    cols = ['avg_depth'] + [f'{fn}_{i+1}' for i in range(n_samples) for fn in ['mean', 'var']]
    coverage_table = pd.DataFrame(coverage, index=contigs, columns=cols)
    coverage_table = pd.concat([ctg_len, coverage_table], axis=1)
    coverage_table.index.name = "contig_id"

    return coverage_table

def metabat2_to_others(table, suffix):

    # Metabat2 format
    table.to_csv(f'coverage_metabat2-{suffix}.tsv', sep="\t")
    # CONCOCT format
    table.iloc[:, 2::2].to_csv(f'coverage_concoct-{suffix}.tsv', sep="\t")
    # Maxbin2 format
    for i, col in enumerate(table.columns[2::2]):
        table[col].to_csv(f"coverage_maxbin2-{suffix}-{i}.tsv", header=True, index=True, sep='\t')
    
def main():
    args = parser()

    if Path(args.abundance).suffix in {'.h5', '.hdf5'}:
        table = h5_to_metabat2(args.abundance)
    else:
        table =  pd.read_csv(args.abundance, sep='\t', index_col=0)

    metabat2_to_others(table, args.suffix)
    
if __name__ == '__main__':
    main()
