#!/usr/bin/env python

import re
import argparse

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt


sns.set(**{'context': 'paper', 'style': 'darkgrid'})
FS = 12

RC = {'axes.labelsize': FS+2,
      'legend.fontsize': FS,
      'axes.titlesize': FS+2,
      'xtick.labelsize': FS-2, 'ytick.labelsize': FS-2}

plt.rcParams.update(**RC)

CATPLOT_KW = dict(sharey='row', kind='box', showfliers=False,
                  margin_titles=True, legend=False, height=2, aspect=1.5)

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=str)
    args = parser.parse_args()

    return args

def load(filename):
    data = pd.read_csv(filename, index_col=0).drop(columns='n_test')
    sim_info = data.index.str.split(r'[\-_]')

    for (name, pos) in [('nb_genomes', 1), ('samples', 3), ('coverage', 5)]:
        data[name] = sim_info.str[pos]
        data[name] = pd.Categorical(
            data[name],
            sorted(data[name].unique(), key=lambda x: int(re.findall(r'\d+', x)[0]))
        )

    main_metrics = ['accuracy', 'F1', 'AUC']
    other_metrics = ['TN', 'TP', 'FN', 'FP']
    meta = ['nb_genomes', 'samples', 'coverage']

    data = dict(
        main=data.drop(columns=other_metrics).melt(id_vars=meta),
        suppl=data.drop(columns=main_metrics).melt(id_vars=meta)
    )

    return data

def plot(data, output=None, ymax=None):

    g = sns.catplot(data=data, x='samples', y='value', hue='coverage',
                    row='variable', col='nb_genomes',
                    **CATPLOT_KW)
    (g.add_legend().despine(left=True)
     .set_titles(col_template='{col_name} genomes', row_template='')
     .set_ylabels(''))

    rownames = data.variable.unique()

    for (row, col, _), facet_data in g.facet_data():
        ax_k = g.facet_axis(row, col)
        if ymax is not None:
            ax_k.set_ylim([ax_k.get_ylim()[0], ymax])

        if col == 0:
            ax_k.set_ylabel(rownames[row])

        if facet_data.empty:
            ax_k.set_axis_off()
            ax_k.set_title('')

    if output is not None:
        plt.savefig(output, bbox_inches='tight')
    else:
        plt.show()

def main():

    args = parse_args()

    data = load(args.results)
    plot(data['main'], output='figures/figure-5_NN-metrics.eps', ymax=1)
    plot(data['suppl'], output='figures/supplementary-nn-metrics.pdf')

if __name__ == '__main__':
    main()
