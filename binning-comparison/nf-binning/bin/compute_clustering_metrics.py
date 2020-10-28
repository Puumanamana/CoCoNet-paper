#!/usr/bin/env python

import argparse
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import SeqIO
import sklearn.metrics


METRICS = ['adjusted_rand_score', 'homogeneity_score', 'completeness_score']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bins', type=str, nargs='+')
    parser.add_argument('--truth', type=str)
    parser.add_argument('--fasta', type=str, default=None)    
    parser.add_argument('--output', type=str)
    parser.add_argument('--name', type=str, default='dataset')
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

    return results.astype(int)


def n_combs(x):
    return x * (x-1) // 2 

def eval_true_bin(assn):
    """
    assn (list): contig assignments for a given true bin
    """

    freqs = assn.value_counts().values
    tp = n_combs(freqs).sum()
    fn = n_combs(len(assn)) - tp
    
    return (tp, fn)

def eval_pred_bin(assn):
    """
    assn (list): contig assignments for a given predicted bin
    """

    freqs = assn.value_counts()
    tp = n_combs(freqs).sum()
    fp = n_combs(len(assn)) - tp

    return fp

def compute_scores(data, output=None, metrics=None, name=None):
    scores = dict()
    
    for method in data.columns.drop('truth'):
        grouped_by_truth = data.groupby('truth')[method].agg(eval_true_bin)
        grouped_by_pred = data.groupby(method).truth.agg(eval_pred_bin)

        scores[method] = dict(
            TP=grouped_by_truth.str[0].sum(),
            FN=grouped_by_truth.str[1].sum(),
            FP=grouped_by_pred.sum(),
        )
        scores[method]['TN'] = n_combs(len(data)) - sum(scores[method].values())

        # Other metrics
        for metric in metrics:
            fn = getattr(sklearn.metrics, metric)
            scores[method][metric] = fn(data.truth, data[method])

    scores = pd.DataFrame(scores).T.astype({col: int for col in ['TP', 'TN', 'FP', 'FN']})
    scores['dataset'] = name

    if output is not None:
        scores.to_csv(output)

    return scores
    

if __name__ == '__main__':
    args = parse_args()
    data = load(args.bins, truth=args.truth)
    print('---- Global scores ----')
    compute_scores(data, output=args.output, metrics=args.metrics, name=args.name)

    if args.fasta is not None:
        scores = {}
        contigs_size = pd.Series({ctg.id: len(ctg.seq)
                                  for ctg in SeqIO.parse(args.fasta, 'fasta')})
        quant_bins = pd.qcut(contigs_size, 5)
        
        for itv in quant_bins.cat.categories:
            subset = data.loc[quant_bins == itv]
            scores[f'[{int(itv.left)}, {int(itv.right)}]'] = compute_scores(subset, output=None, metrics=args.metrics)

        scores = pd.concat(scores).reset_index(level=1, drop=True)

        print(scores)
