#!/usr/bin/env python3

from pathlib import Path
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str)
    parser.add_argument('--bins', type=str)
    args = parser.parse_args()

    return args

def merge_contigs(fasta, assignments):
    contigs = {bin_id: '' for bin_id in set(assignments.values())}
    bin_max = max(contigs.keys())

    with open(fasta) as handle:
        for (ctg_descr, seq) in SimpleFastaParser(handle):
            ctg_id = ctg_descr.split()[0]
            if ctg_id in assignments:
                contigs[assignments[ctg_id]] += seq
            else:
                # singleton bin
                bin_max += 1
                contigs[bin_max] = seq
                assignments[ctg_id] = bin_max

    return contigs

def write_bins_fa(contigs, output):
    with open(output, 'w') as handle:
        for (bin_id, seq) in sorted(contigs.items(), key=lambda x: int(x[0])):
            handle.write(f'>bin_{bin_id}\n{seq}\n')

def write_bins_csv(assignments, output):
    with open(output, 'w') as handle:
        for (ctg, bin_id) in assignments.items():
            handle.write(f'{ctg},{bin_id}\n')

def main():
    '''
    '''

    args = parse_args()
    assignments = pd.read_csv(args.bins, header=None, names=['contig', 'bin_id'])
    assignments.bin_id = pd.factorize(assignments.bin_id)[0]
    assignments = assignments.set_index('contig').bin_id.to_dict()

    prefix = Path(args.bins).stem
    merged = merge_contigs(args.fasta, assignments)
    write_bins_fa(merged, f'{prefix}-merged.fasta')
    write_bins_csv(assignments, f'{prefix}-complete.csv')

if __name__ == '__main__':
    main()

