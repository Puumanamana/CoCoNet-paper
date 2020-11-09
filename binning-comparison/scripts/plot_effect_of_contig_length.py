#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import sklearn.metrics
import pingouin as pg

import matplotlib.pyplot as plt
import seaborn as sns


sns.set(**{'context': 'paper', 'style': 'darkgrid'})
FS = 12

RC = {'axes.labelsize': FS+2,
      'legend.fontsize': FS,
      'axes.titlesize': FS+2,
      'xtick.labelsize': FS-2, 'ytick.labelsize': FS-2}

plt.rcParams.update(**RC)

METRICS = ['adjusted_rand_score', 'homogeneity_score', 'completeness_score']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root-dir', type=str)
    args = parser.parse_args()

    return args

def find_files(folder):
    mapping = []
    for bin_file in Path(folder, 'binning/CoCoNet').glob('coconet-*.csv'):
        fasta = Path(bin_file.parent, bin_file.stem, 'assembly-filtered.fasta')
        info = re.split(r'[_\-]', Path(bin_file).stem)
        mapping.append([info[2], info[4], info[6], info[7], bin_file, fasta])

    return pd.DataFrame(mapping, columns=['genomes', 'samples', 'coverage', 'replicate', 'bins', 'fasta'])
        

def compute_scores(bin_file, fasta, intervals):
    data = pd.read_csv(bin_file, index_col=0, header=None, names=['contig', 'pred'])
    data['truth'] = pd.factorize(data.index.str.split('|').str[0])[0]

    contigs_size = pd.Series({ctg.id: len(ctg.seq) for ctg in SeqIO.parse(fasta, 'fasta')})
    size_bins = pd.cut(contigs_size, bins=intervals)

    scores = {}
    for itv in size_bins.cat.categories:
        subset = data.loc[size_bins == itv]
        itv_str = f'{int(itv.left)}-{int(itv.right)}'
        scores[itv_str] = [getattr(sklearn.metrics, metric)(subset.truth, subset.pred)
                           for metric in METRICS]

    scores = pd.DataFrame(scores, index=METRICS)

    return scores
    
if __name__ == '__main__':
    args = parse_args()

    files_info = find_files(args.root_dir)
    intervals = [2000, 3000, 10000, 100000]

    scores = []
    for _, entry in files_info.iterrows():
        score = (compute_scores(entry.bins, entry.fasta, intervals)
                 .assign(samples=entry.samples, coverage=entry.coverage,
                         genomes=entry.genomes, replicate=entry.replicate))
        scores.append(score)

    scores = pd.concat(scores).set_index(['samples', 'coverage', 'genomes', 'replicate'],
                                         append=True)
    scores.columns.name = 'interval'
    scores = scores.stack().loc['adjusted_rand_score'].rename('ARI').reset_index()

    print('====== Kruskal-Wallis test ======')
    print(pg.kruskal(data=scores, dv='ARI', between='interval'))

    g = sns.catplot(data=scores, x='interval', y='ARI', kind='box', showfliers=False)
    g.set(xlabel='Contig length category')

    plt.savefig('figures/supplementary-effect-of-contig-length.pdf', bbox_inches='tight')

