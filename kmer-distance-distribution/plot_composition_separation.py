import tempfile
import argparse
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from sklearn.metrics.pairwise import paired_distances
import seaborn as sns
import matplotlib.pyplot as plt

from coconet.fragmentation import make_pairs
from coconet.core.generators import CompositionGenerator

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', type=str, default='/home/cedric/db/viral.genomic.ACGT.fasta')
    parser.add_argument('--fragment-lengths', type=int, nargs='+', default=[256, 512, 1024, 2048])
    parser.add_argument('--k', type=int, default=4)
    args = parser.parse_args()

    args.db = Path(args.db)

    return args

def compute_distances(pair_file, fasta, kmer=4, rc=True, norm=True, batch_size=None):
    kmer_generator = CompositionGenerator(
        pair_file, fasta, batch_size=batch_size,
        kmer=kmer, rc=rc, norm=norm)

    x_compo = map(lambda x: x.numpy(), next(kmer_generator))
    distances = paired_distances(*x_compo, 'cosine')

    return distances

def run(fl, k, db, db_path):
    print('Processing fragment length = {}'.format(fl))
    pairs_out = tempfile.mkdtemp() + '/pairs.npy'

    pairs = make_pairs(db, 64, fl, n_examples=len(db)*5, output=pairs_out)
    truth = ['Within species' if sp1 == sp2 else 'Across species' for sp1, sp2 in pairs['sp']]
    distances = pd.DataFrame({
        'dist': compute_distances(pairs_out, db_path, kmer=k, batch_size=len(pairs)),
        'truth': truth,
        'fl': fl})

    shutil.rmtree(Path(pairs_out).parent)
    return distances

def plot(distances):

    # sns.boxplot(x="fl", y="dist", hue="truth", data=distances, showfliers=False, hue_order=['Across species', 'Within species'])
    sns.violinplot(x="fl", y="dist", hue="truth", data=distances,
                   hue_order=['Across species', 'Within species'],
                   split=True, inner=None, gridsize=500)
    sns.despine(left=True)

    # plt.ylim([0, 0.2])
    plt.xlabel('Fragment length (bp)', fontsize=14)
    plt.ylabel('Distance', fontsize=14)
    plt.legend(title='')

    plt.savefig('Figure 1-Composition_separation_with_fl.pdf', transparent=True)
    plt.show()


def main():
    args = parse_args()

    db_genomes = [(seq.id, str(seq.seq)) for seq in SeqIO.parse(args.db, 'fasta')
                  if len(seq.seq) >= max(args.fragment_lengths)]
    print('{} genomes'.format(len(db_genomes)))

    distances_all = []

    results = [run(fl, args.k, db_genomes, args.db) for fl in args.fragment_lengths]
    distances_all = pd.concat(results)

    plot(distances_all)

if __name__ == '__main__':
    main()
