import argparse
from pathlib import Path
import re

import pandas as pd
import numpy as np
import h5py


def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('root_dir', type=str)
    parser.add_argument('--filtering', action='store_true')
    parser.add_argument('--summary', action='store_true')
    args = parser.parse_args()

    return args

def parse_sim_name(name):
    info = re.split(r'[_\-]', name)
    # (genomes, samples, coverage, replicate)
    return [info[1], info[5][:-1], info[3]]

def sim_stats(root, sim, noise=0.1):
    h5_file = Path(root, 'preprocessing', sim, 'coverage.h5')
    with h5py.File(h5_file, 'r') as handle:
        prevalence = np.mean([sum(cov[:].mean(axis=1) > noise) for cov in handle.values()])
        n_contigs = len(handle)
        bins = pd.factorize([x.split('|')[0] for x in handle.keys()])[0]

    infos = parse_sim_name(sim) + [n_contigs, np.bincount(bins).mean(), prevalence]

    return infos

def sim_filtering_summary(root, sim):

    # log file for raw number of contigs
    log_file = Path(root, 'preprocessing', sim, 'coconet.log')
    n_ctg = int(re.findall('[0-9,]+', next(open(log_file)))[-1].replace(',', ''))

    # exclude.tsv for length filtered contigs
    length_file = Path(log_file.parent, 'exclude.tsv')
    length_filtered = {x.split('\t')[0] for i, x in
                       enumerate(open(length_file)) if i > 0}

    # other exclude.tsv for prevalence filtered contigs
    prev_file = Path(root, 'binning', 'CoCoNet', f'coconet-{sim}', 'exclude.tsv')

    if not prev_file.is_file():
        return

    prev_filtered = {x.split('\t')[0] for i, x in
                     enumerate(open(prev_file)) if i > 0}

    # dtr.tsv for complete genomes filtering
    dtr_file = Path(prev_file.parent, 'dtr.tsv')
    dtr_filtered = {x.split('\t')[0] for x in open(dtr_file)}

    filters = [n_ctg, n_ctg - len(length_filtered),
               n_ctg - len(length_filtered) - len(prev_filtered.difference(dtr_filtered))]

    return parse_sim_name(sim) + filters

def main():

    args = parse_args()

    summary_list = []
    for i, sim in enumerate(Path(args.root_dir, 'preprocessing').glob('*')):
        if args.filtering:
            entry = sim_filtering_summary(args.root_dir, sim.name)
        elif args.summary:
            print(f'Processed {i}/80', end='\r')
            entry = sim_stats(args.root_dir, sim.name)

        if entry is not None:
            summary_list.append(entry)

    index_cols = ['Number of genomes', 'Coverage', 'Number of samples']
    if args.filtering:
        value_cols = ['raw', 'L>2048bp', 'prevalence>1']
    else:
        value_cols = ['Number of contigs', 'Average bin size', 'Average prevalence']

    summaries = pd.DataFrame(summary_list, columns=index_cols + value_cols)

    for idx in index_cols:
        suffix = ''
        if 'coverage' in idx.lower():
            suffix = 'X'

        summaries[idx] = pd.Categorical(
            summaries[idx] + suffix,
            [f'{x}{suffix}' for x in sorted(summaries[idx].unique(), key=int)]
        )

    if args.filtering:
        summaries = (summaries
                     .groupby(index_cols)
                     .mean()
                     .applymap(lambda x: f'{x:,.0f}')
                     .reset_index())
        summaries.index = [f'Sim-{i+1}' for i in range(len(summaries))]
        summaries.index.name = 'Dataset'
        summaries.to_csv('simulations-filtering.tsv', sep='\t')

    else:
        summaries = (summaries
                     .groupby(index_cols)
                     .agg(['mean', 'std']).stack(level=0))
        summaries['txt'] = [f'{mu:.2f} +-{sd:.2f}' for (mu, sd) in summaries.values]

        summaries = summaries['txt'].unstack().reset_index()
        summaries.index = [f'Sim-{i+1}' for i in range(len(summaries))]
        summaries.index.name = 'Dataset'

        summaries.to_csv('simulations-stats.tsv', sep='\t')

if __name__ == '__main__':
    main()
