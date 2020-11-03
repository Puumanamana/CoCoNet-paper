import re
import sys

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

data = pd.read_csv(sys.argv[1]).set_index(['genomes', 'samples', 'coverage',
                                           'replicate', 'mode', 'last_batch'])
data = (data
        .drop(columns=['TP', 'TN', 'FP', 'FN'])
        .rename_axis(columns='metric').stack().rename('score')
        .reset_index())

data['setting'] = pd.Categorical(data.samples.astype(str) + ' samples ('
                                 + data.coverage.astype(str) + ')')
data.setting = data.setting.cat.reorder_categories(
    sorted(data.setting.cat.categories, key=lambda x: [int(x) for x in re.findall(r'\d+', x)])
)

data = data[data.metric == 'AUC'].rename(columns={'score': 'AUC'}).replace(dict(with_var='variable_coverage', no_var='flat_coverage'))

grid = sns.catplot(
    kind='box', data=data, x='mode', y='AUC', showfliers=False,
    col='coverage', row='samples', sharey=False, margin_titles=True,
)

grid.set_titles(row_template='{row_name} samples')
grid.savefig('performance-with-coverage-variability.pdf')

plt.show()
