#!/usr/bin/env python

import re
import argparse

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt


sns.set(**{'context': 'paper', 'style': 'darkgrid'})
FS = 10

RC = {'axes.labelsize': FS+2,
      'legend.fontsize': FS,
      'axes.titlesize': FS+2,
      'xtick.labelsize': FS-2, 'ytick.labelsize': FS-2}

plt.rcParams.update(**RC)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=str)
    parser.add_argument('--sim', action='store_true')    
    args = parser.parse_args()

    return args

def parse_sim(data):
    sim_info = data.index.str.split(r'[\-_]')

    data['nb_genomes'] = sim_info.str[2].astype(int)
    data['replicate'] = sim_info.str[7].astype(int)
    data['setting'] = pd.Categorical(sim_info.str[4] + ' samples (' + sim_info.str[6] + ')')

    data.setting = data.setting.cat.reorder_categories(
        sorted(data.setting.cat.categories, key=lambda x: [int(x) for x in re.findall('\d+', x)])
    )
    data = data.reset_index(drop=True)

    data = data.melt(id_vars=['nb_genomes', 'setting', 'method', 'replicate'],
                     var_name='metric', value_name='score')

    return data
    

def plot_sim(data, output=None):
    graph = sns.catplot(data=data, kind='bar', sharey=False, margin_titles=True, legend=False,
                        x='metric', y='score', hue='method',
                        order=['adjusted_rand_score', 'completeness_score','homogeneity_score'],
                        col='nb_genomes', row='setting')
    graph.set_titles(row_template='{row_name}', col_template='{col_name} genomes')
    graph.set(xlabel='', ylabel='')
    graph.add_legend(title='')

    if output is None:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')

def plot_SA(data, output=None):
    
    graph = sns.catplot(data=data, kind='bar', sharey=False, margin_titles=True, legend=False,
                        x='metric', y='score', hue='method',
                        order=['adjusted_rand_score', 'completeness_score','homogeneity_score'])
    # graph.set_titles(row_template='{row_name}', col_template='{col_name} genomes')
    graph.set(xlabel='', ylabel='')
    graph.add_legend(title='')

    if output is None:
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')

def main():

    args = parse_args()

    data = pd.read_csv(args.results, index_col=0).drop(columns='dataset')
    data['method'] = data.index.str.split(r'-|_').str[0]
    data = data.replace(dict(coconet='CoCoNet', metabat2='Metabat2', concoct='CONCOCT'))
    data['method'] = pd.Categorical(data.method, ['CONCOCT', 'Metabat2', 'CoCoNet']) 

    if args.sim:
        data = parse_sim(data)
        plot_sim(data, output='figure-6.pdf')

    else:
        data = (data.drop(columns=['FN', 'FP', 'TN', 'TP'])
                .melt(var_name='metric', value_name='score', id_vars='method'))
        plot_SA(data, output='figure-7A.pdf')

if __name__ == '__main__':
    main()