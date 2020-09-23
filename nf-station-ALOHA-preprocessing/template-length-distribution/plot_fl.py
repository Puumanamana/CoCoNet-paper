import os
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_last_line(f):
    with open(f, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        last_line = f.readline().decode()

    return last_line

def get_freqs(bins=20):
    data = {}
    for filename in Path('.').glob('tlen-*.txt'):
        n = int(get_last_line(filename).strip().split()[-1])
        name = filename.stem.split('-')[-1]
    
        data[name] = np.zeros(n+1)
        for entry in open(filename):
            freq, tlen = entry.split()
            data[name][int(tlen)] = int(freq)

    for f, counts in data.items():
        counts = pd.Series(counts[0:2000]).reset_index()
        counts.columns = ['tlen', 'freq']
        counts['sample'] = f
        
        data[f] = counts

    data = pd.concat(data.values())
    data.tlen = pd.qcut(data.tlen, q=bins)
    data.tlen = data.tlen.cat.rename_categories({x: f'{int(x.left)}-{int(x.right)}' for x in data.tlen.values})
    
    return data

def plot_tlen(size=1, margin=0.2):
    counts = get_freqs()

    grid = sns.FacetGrid(data=counts, col='sample', sharey=False, col_order=sorted(counts['sample'].unique()), height=size)
    grid.map(sns.barplot, 'tlen', 'freq', order=counts.tlen.drop_duplicates())
    grid.set_xticklabels(rotation=90)
    plt.subplots_adjust(bottom=margin)
    grid.savefig('template_length.pdf')

if __name__ == '__main__':
    plot_tlen(2, 0.1)
    import ipdb;ipdb.set_trace()
