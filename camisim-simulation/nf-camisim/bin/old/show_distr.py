import pandas as pd
import numpy as np

from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt

seqd = 1e9*0.05762946

def get_lg_params(x):
    shape, loc, scale = stats.lognorm.fit(x, loc=0)
    return (np.log(scale), shape)

grieg = pd.read_csv('Grieg.csv', sep=' ', header=None)*seqd
grieg.columns = [f'g{i}' for i in range(grieg.shape[1])]

all_distr = pd.concat([pd.read_csv(f'distribution_{i}.txt', sep='\t',
                                   header=None, index_col=0, names=['ctg', f'd{i}'])
                       for i in range(4)], axis=1)
all_distr *= seqd

data = pd.concat([grieg.melt(), all_distr.melt()])
data = data[data.value>0]

prms = data.groupby('variable').agg(get_lg_params)
print(prms)

data.value = np.log(data.value)


g = sns.FacetGrid(data, hue='variable')
g.map(sns.kdeplot, 'value').add_legend()
plt.show()
