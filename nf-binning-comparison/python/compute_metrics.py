import argparse
from pathlib import Path

import pandas as pd
import sklearn.metrics

METRICS = ['adjusted_rand_score', 'homogeneity_score', 'completeness_score']

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--bins', type=str, nargs='+')
    parser.add_argument('--truth', type=str)
    parser.add_argument('--output', type=str)
    args = parser.parse_args()

    return args


def compute_score(metric, bin_file, truth_file=None):
    bins = pd.read_csv(bin_file, header=None, names=['contig', 'bin'])

    if truth_file is None: # e.g. truth is in the contig name 
        truth = pd.factorize(bins.index.str.split('_').str[0])[0]
    else: # real dataset
        truth = pd.read_csv(truth_file, header=None, names=['contig', 'bin'])
    
    score = getattr(sklearn.metrics, metric)(truth, bins)

    return score


def compute_scores_all_metrics():
    '''
    '''

    args = parse_args()

    (ds, tool) = Path(args.bins).stem.split('-')

    results = []
    for metric in METRICS:
        score = compute_score(metric, args.bins, args.truth)
        results.append([ds, tool, metric, score])

    results = pd.DataFrame(results, columns=['dataset', 'tool', 'metric', 'score'])
    

if __name__ == '__main__':
    compute_scores_all_metrics()
