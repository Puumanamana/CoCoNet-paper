#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd
import sklearn.metrics


METRICS = ['accuracy_score', 'roc_auc_score', 'f1_score']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', type=str)
    parser.add_argument('--name', type=str)    
    parser.add_argument('--output', type=str)
    args = parser.parse_args()

    args.test = Path(args.test)
    
    return args

def compute_scores(truth, pred):
    pred_bin = (pred > 0.5).astype(int)
    
    matrix = sklearn.metrics.confusion_matrix(truth, pred_bin)

    scores = pd.Series(dict(
        n_test=matrix.sum(),
        TN=matrix[0, 0],
        TP=matrix[1, 1],
        FP=matrix[0, 1],
        FN=matrix[1, 0],
        accuracy=sklearn.metrics.accuracy_score(truth, pred_bin),
        F1=sklearn.metrics.f1_score(truth, pred_bin),
        AUC=sklearn.metrics.roc_auc_score(truth, pred)
    ))

    return scores
    
def main():
    args = parse_args()
    
    results = pd.read_csv(args.test)
    scores = compute_scores(results.truth, results.combined)
    scores.name = args.name

    scores = pd.DataFrame(scores).T
        
    scores = scores.astype(
        dict([(col, float) if col in {'AUC', 'F1', 'accuracy'} else (col, int)
              for col in scores.columns])
    )

    if args.output is None:
        print(scores)

    else:
        scores.to_csv(args.output)
        
    
if __name__ == '__main__':
    main()

