#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO


sns.set(**{'context': 'paper', 'style': 'darkgrid'})
FS = 22

RC = {'axes.labelsize': FS+2,
      'legend.fontsize': FS,
      'axes.titlesize': FS+2,
      'xtick.labelsize': FS-2, 'ytick.labelsize': FS-2}

plt.rcParams.update(**RC)

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=str)
    parser.add_argument('--refseq', type=str, help='RefSeq database in fasta format')
    args = parser.parse_args()

    return args

def load(folder):
    
    summary_files = Path(folder).glob('*/quality_summary.tsv')

    summaries = {}

    for filename in summary_files:
        name = filename.parent.name.split('-')[-1]
        
        summary = pd.read_csv(
            filename, index_col=0, sep='\t',
            usecols=['contig_id', 'contig_length', 'viral_genes', 'checkv_quality',
                     'miuvig_quality', 'completeness', 'contamination'],
            low_memory=False
        ).rename(columns={'contig_length': 'bin_size'})
        summaries[name] = summary

    summaries = pd.concat(summaries)
    summaries.index.names = ['method', 'bin_id']
    grouped = summaries.groupby(level=0)

    return dict(
        summary=summaries,
        checkv_quality=grouped.checkv_quality.value_counts(),
        bin_size=grouped.bin_size.mean(),
        miuvig_quality=grouped.miuvig_quality.value_counts(),
        mean_contamination=grouped.contamination.mean().apply(lambda x: f'{x:.2f}%')
    )

def main():

    args = parse_args()


    assembly_quality = load(args.folder)

    plot_data = (
        assembly_quality['checkv_quality'].unstack()
        .assign(bin_size=0)
        .stack()
        .rename('count').reset_index()
    )

    plot_data.method = pd.Categorical(plot_data.method.replace('metabat2', 'Metabat2'),
                                      ['CONCOCT', 'Metabat2', 'CoCoNet'])
    grid = sns.catplot(data=plot_data, kind='bar', x='method', y='count', edgecolor='black',
                       col='checkv_quality', sharey=False, sharex=False, col_wrap=2, aspect=2)
    grid.set_titles(col_template='#{col_name} bins')
    grid.set(xlabel='', ylabel='Count')

    for ax in [1, 3]:
        grid.facet_axis(0, ax).set_ylabel('')

    last_ax = grid.facet_axis(0, 5)

    bin_sizes = assembly_quality['summary'].bin_size.reset_index()
    bin_sizes.method = pd.Categorical(bin_sizes.method.replace('metabat2', 'Metabat2'),
                                      ['CONCOCT', 'Metabat2', 'CoCoNet'])

    if args.refseq is not None:
        ctg_size_refseq = pd.DataFrame([
            ['RefSeq', seq.id, len(seq.seq)]
            for seq in SeqIO.parse(args.refseq, 'fasta')
        ], columns=['method', 'bin_id', 'bin_size'])
        bin_sizes = pd.concat([bin_sizes, ctg_size_refseq])
        bin_sizes.method = pd.Categorical(
            bin_sizes.method, ['CONCOCT', 'Metabat2', 'CoCoNet', 'RefSeq']
        )

    sns.boxplot(data=bin_sizes, x='method', y='bin_size', ax=last_ax, showfliers=False)
    last_ax.set_yscale('log')
    last_ax.set_title('Distribution of bin sizes')
    last_ax.set_xlabel('')
    last_ax.set_ylabel('Size (bp)')

    plt.savefig('figures/figure-7B-checkv-SA.eps', bbox_inches='tight')

if __name__ == '__main__':
    main()
