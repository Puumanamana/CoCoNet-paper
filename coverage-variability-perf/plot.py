from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = []
for csv in Path('results').glob('*.csv'):
    (n_genomes, coverage, n_samples, iteration) = csv.stem.split('-')[-1].split('_')

    data_i = pd.read_csv(csv).assign(
        n_genomes=n_genomes,
        coverage=coverage,
        n_samples=n_samples,
        iteration=iteration
    )

    data.append(data_i)

data = pd.concat(data).reset_index(drop=True)

data_last = data[data.batch==3800]
data_plot = data_last.melt(
    id_vars=[x for x in data.columns if x not in {'acc', 'AUC', 'prec', 'recall'}]
)

sns.catplot(kind='box', data=data_plot, x='variable', y='value', hue='mode', col='coverage', row='n_samples')
plt.show()
