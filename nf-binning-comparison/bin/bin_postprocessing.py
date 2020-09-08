#!/usr/bin/env python3

import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser

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
        for (ctg_descr, ctg) in SimpleFastaParser(handle):
            ctg_id = ctg_descr.split()[0]
            if ctg_id in assignments:
                contigs[assignments[ctg_id]] += ctg
            else:
                # singleton bin
                bin_max += 1
                contigs[bin_max] = ctg
                assignments[ctg_id] = bin_max

    return contigs

def write_bins_fa(contigs):
    with open("merged.fasta", 'w') as handle:
        for (bin_id, seq) in sorted(contigs.items(), key=lambda x: int(x[0])):
            handle.write(f'>bin_{bin_id}\n{seq}\n')

def write_bins_csv(assignments):
    with open('bins_with_singletons.csv', 'w') as handle:
        for (ctg, bin_id) in assignments.items():
            handle.write(f'{ctg},{bin_id}\n')

def main():
    '''
    '''

    args = parse_args()
    assignments = {ctg: int(b) for (ctg, b) in map(lambda x: x.strip().split(','), open(args.bins))}
    
    merged = merge_contigs(args.fasta, assignments)
    write_bins_fa(merged)
    write_bins_csv(assignments)

if __name__ == '__main__':
    main()

