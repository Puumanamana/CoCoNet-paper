#!/usr/bin/env python3

import argparse

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


sns.set(**{'context': 'paper', 'style': 'darkgrid'})
FS = 14

RC = {'axes.labelsize': FS+2,
      'legend.fontsize': FS,
      'axes.titlesize': FS+2,
      'xtick.labelsize': FS-2, 'ytick.labelsize': FS-2}

plt.rcParams.update(**RC)

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('table')
    args = parser.parse_args()

    return args

def main():
    '''
    '''

    args = parse_args()

    data = pd.read_csv(args.table, na_values='-')
    data['process'] = data.process.str.split(':').str[-1]

    processes = {'COCONET_RUN', 'CONCOCT', 'METABAT2'}
    data = data.loc[data.process.isin(processes)]

    data['n_contigs'] = data.tag.str.split('-').str[2].astype(int)
    data['time'] = pd.to_timedelta(data.realtime).dt.seconds / 60

    rss_values = data.peak_rss.str.split(' ', expand=True)
    data['rss'] = rss_values[0].astype(float) / rss_values[1].replace(dict(MB=2**10, GB=1))

    data = (data
            .replace(dict(COCONET_RUN='CoCoNet', METABAT2='Metabat2'))
            .set_index(['n_contigs', 'process'])[['rss', 'time']]
            .rename_axis(columns='metric')
            .stack().rename('score').reset_index())

    g = sns.relplot(data=data, col='metric', hue='process', style='process',
                    facet_kws=dict(sharey=False),
                    err_style='bars', ci='sd', linewidth=2,
                    x='n_contigs', y='score', kind='line', markers=True, dashes=False)

    g.set(xlabel='Number of contigs')
    g.axes[0, 0].set_title('Memory Usage')
    g.axes[0, 0].set_ylabel('Peak memory (GB)')
    g.axes[0, 1].set_title('Running time')
    g.axes[0, 1].set_ylabel('Time (min)')
    g._legend.set_title(None)

    g.savefig('time-memory.pdf')
    plt.show()

if __name__ == '__main__':
    main()
