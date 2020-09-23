#!/usr/bin/env python

import argparse
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import sklearn.metrics


METRICS = ['adjusted_rand_score', 'homogeneity_score', 'completeness_score']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bins', type=str, nargs='+')
    parser.add_argument('--truth', type=str)
    parser.add_argument('--output', type=str)
    parser.add_argument('--metrics', type=str, nargs='+', default=METRICS)    
    args = parser.parse_args()

    return args


def load_bins(filename):
    if not Path(filename).is_file():
        raise FileNotFoundError(f"{filename} not found")

    df = pd.read_csv(filename, index_col=0, header=None,
                     names=['contigs', 'clusters'])
    df.clusters = pd.factorize(df.clusters)[0]
        
    return df


def load(files, truth=None):
    results = (pd.concat({Path(f).stem: load_bins(f) for f in files})
               .unstack(level=0)
               .clusters)

    if truth is None or not Path(truth).is_file():
        # Truth in contig name
        results['truth'] = pd.factorize(results.index.str.split('|').str[0])[0]
    else:
        # Load truth file
        truth = pd.read_csv(truth, index_col=0, header=None,
                            names=['contigs', 'clusters'])
        truth.clusters = pd.factorize(truth.clusters)[0]
        results['truth'] = truth.reindex(results.index)
        results.dropna(how='any', inplace=True)

    if results.size == 0:
        raise RuntimeError('No contigs shared between truth and all methods')

    return results

def compute_scores(data, output=None, metrics=None):
    scores = {}
    for method in data.columns.drop('truth'):
        scores[method] = {}
        for metric in metrics:
            fn = getattr(sklearn.metrics, metric)
            scores[method][metric] = fn(data.truth, data[method])

        truth_npy = data.truth.values
        pred_npy = data[method].values
        matrix = np.zeros((2, 2))

        for (i, j) in combinations(range(len(truth_npy)), 2):
            matrix[
                int(truth_npy[i] == truth_npy[j]),
                int(pred_npy[i] == pred_npy[j])
            ] += 1
        
        scores[method]['TN'] = matrix[0, 0]
        scores[method]['TP'] = matrix[1, 1]        
        scores[method]['FP'] = matrix[0, 1]
        scores[method]['FN'] = matrix[1, 0]        

    scores = pd.DataFrame(scores).T

    if output is None:
        print(scores)
    else:
        scores.to_csv(output)
    

if __name__ == '__main__':
    args = parse_args()
    data = load(args.bins, truth=args.truth)
    compute_scores(data, output=args.output, metrics=args.metrics)
