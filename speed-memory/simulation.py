#!/usr/bin/env python3

import argparse
import random
from pathlib import Path

import numpy as np
from Bio import SeqIO

import h5py


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--db', type=str)
    parser.add_argument('--n-contigs', type=int, default=10)
    parser.add_argument('--n-samples', type=int, default=5)
    parser.add_argument('--min-genome-len', type=int, default=4000)
    parser.add_argument('--replicate', type=int, default=1)    
    args = parser.parse_args()

    return args

def split_genomes(genomes, coverage, n_contigs=10, suffix='',
                  min_ctg_len=2048):

    outdir = Path('simulations', f'sim-{suffix}')
    outdir.mkdir(exist_ok=True, parents=True)

    ctg_handle = open(f'{outdir}/assembly.fasta', 'w')
    h5_handle = h5py.File(f'{outdir}/coverage_contigs.h5', 'w')

    genome_count = 0
    ctg_count = 0

    while True:
        genome = random.choice(genomes)
        genome_seq = str(genome.seq)
        n_splits = 1+random.choice(range(len(genome.seq)//min_ctg_len))

        starts = sorted(np.random.randint(0, len(genome_seq)-min_ctg_len, n_splits))
        ends = np.diff(starts, append=len(genome.seq))
        itv = filter(lambda x: (x[1]-x[0])>=2048, zip(starts, ends))

        for i, (start, end) in enumerate(itv):
            ctg_handle.write(f'>V{genome_count}|{i}\n{genome_seq[start:end]}\n')
            h5_handle.create_dataset(
                f'V{genome_count}|{i}',
                data=coverage[genome.id][:, start:end]
            )
            ctg_count += 1

            if ctg_count >= n_contigs:
                ctg_handle.close()
                h5_handle.close()
                return

        genome_count += 1

def generate_coverage(genomes, n_samples=5):
    coverage = {}
    for genome in genomes:
        mus = np.random.lognormal(1, 2, size=n_samples)
        coverage[genome.id] = np.stack([
            np.random.poisson(mu, len(genome.seq))
            for mu in mus
        ])
    return coverage

def main():
    args = parse_args()

    label = f'{args.n_samples}-{args.n_contigs}-{args.replicate}'

    genomes = [genome for genome in SeqIO.parse(args.db,'fasta')
               if len(genome.seq) > args.min_genome_len]

    coverage = generate_coverage(genomes, n_samples=args.n_samples)

    split_genomes(genomes, coverage, suffix=label, n_contigs=args.n_contigs, min_ctg_len=2048)

if __name__ == '__main__':
    main()
