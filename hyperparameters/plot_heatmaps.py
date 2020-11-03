from itertools import combinations
import argparse

import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


sns.set(**{'context': 'paper', 'style': 'darkgrid'})
FS = 12

RC = {'axes.labelsize': FS,
      'figure.titlesize': FS+20,
      'legend.fontsize': FS,
      'axes.titlesize': FS,
      'xtick.labelsize': FS-2, 'ytick.labelsize': FS-2}

plt.rcParams.update(**RC)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str)
    parser.add_argument('--rows', type=str)
    parser.add_argument('--cols', type=str)
    parser.add_argument('--values', type=str,
                        choices=['completeness', 'homogeneity', 'ARI'], nargs='+')
    parser.add_argument('--plot-all', action='store_true')
    args = parser.parse_args()

    return args

def facetgrid_heatmap(data, rows, cols, values):
    grouped = (data.groupby([rows, cols])[values].mean()
               .rename_axis(columns='metrics')
               .stack().rename('score')
               .reset_index())

    def facet(data, **kwargs):
        data = data.pivot(index=rows, columns=cols, values='score')
        g = sns.heatmap(data, **kwargs)

    g = sns.FacetGrid(grouped, col='metrics')
    cbar_ax = g.fig.add_axes([.92, .3, .02, .4])

    g = g.map_dataframe(facet, cbar_ax=cbar_ax, cmap='hot',
                        linewidths=.1, linecolor='k')
    g.set_titles(col_template="{col_name}", fontweight='bold')
    g.set(xlabel=cols, ylabel=rows)

    for ax in g.axes.flat:
        ax.set_xlabel(ax.get_xlabel().capitalize().replace('_', ' '))
        ax.set_ylabel(ax.get_ylabel().capitalize().replace('_', ' '))

        title = ax.get_title()
        if title != 'ARI':
            ax.set_title(title.capitalize(), fontweight='bold')

    g.fig.subplots_adjust(right=.9, wspace=0.2, bottom=0.2, left=0.2)

    plt.savefig(f'figures/heatmap-{rows}-vs-{cols}.pdf')

def run(csv, *args):

    data = pd.read_csv(csv).drop(columns=['Trial-ID', 'Status', 'Iteration']).round(4)
    data = data.rename(columns=dict(Objective='ARI'))
    facetgrid_heatmap(data, *args)

if __name__ == '__main__':
    args = parse_args()

    if not args.plot_all:
        run(args.data, args.rows, args.cols, args.values)

    else:
        targets = ['ARI', 'homogeneity', 'completeness']
        factors = ['gamma1', 'gamma2', 'theta', 'max_neighbors']

        for (k1, k2) in combinations(factors, 2):
            run(args.data, k1, k2, targets)
