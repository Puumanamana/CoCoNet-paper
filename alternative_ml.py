import argparse
from pathlib import Path
import subprocess

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier

from coconet.dl_util import load_data, get_labels

root = 'paper_data/output_data'

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--sim', type=str, default='2000_3_4_0')
    parser.add_argument('--load-batch', type=int, default=100)
    args = parser.parse_args()

    sim_dir = Path(root, args.sim)
    args.pairs = dict(train=Path(sim_dir, 'pairs_train.npy'),
                      test=Path(sim_dir, 'pairs_test.npy'))
    args.assembly = Path(sim_dir, 'assembly_filtered.fasta')
    args.coverage = Path(sim_dir, 'coverage_filtered.h5')
    args.coconet_result = Path(sim_dir, 'CoCoNet_test.csv')

    if not args.coverage.is_file():
        gz_cov = Path(sim_dir, 'coverage_filtered.h5.gz')
        if not gz_cov.is_file():
            print('Aborting. No hdf5 coverage file found')
            exit(1)
        subprocess.check_output(['unpigz', '-p', '10', str(args.coverage)])
            

    return args

def extract_coconet_score(filename):
    results = pd.read_csv(filename)

    results.combined = (results.combined > 0.5).astype(int)
    score = (results.combined == results.truth).mean()

    return score

def ml_predict(assembly, coverage, pairs, load_batch=1000):
    '''
    '''

    (x_test, x_train_gen) = (load_data(assembly, coverage, pairs[mode], mode=mode,
                                       kmer=4, rc=True, norm=False, wsize=64, wstep=32,
                                       load_batch=load_batch, batch_size=load_batch)
                             for mode in ['test', 'train'])

    (xt11, xt12, xt21, xt22) = [x.detach().numpy() for x in x_test[0] + x_test[1]]
    x_test = np.hstack([xt11, xt12, xt21.reshape(len(xt11), -1), xt22.reshape(len(xt11), -1)])
    
    (y_train, y_test) = (get_labels(pairs['train']),
                         get_labels(pairs['test']).detach().numpy().astype(int)[:, 0])

    clf = RandomForestClassifier(n_jobs=10, warm_start=True, n_estimators=1)

    for i, (X1, X2) in enumerate(x_train_gen, 1):
        y_b = y_train[(i-1)*load_batch:i*load_batch, 0]
        
        (x11, x12, x21, x22) = [x.detach().numpy() for x in X1+X2]

        x_b = np.hstack([x11, x12,
                         x21.reshape(load_batch, -1), x22.reshape(load_batch, -1)])

        clf.fit(x_b,y_b)
        clf.n_estimators += 1

    y_pred = np.argmax(clf.predict_proba(x_test), axis=1)

    return np.mean(y_pred == y_test)


def main():
    args = parse_args()

    ml_score = ml_predict(args.assembly, args.coverage, args.pairs, args.load_batch)
    coco_score = extract_coconet_score(args.coconet_result)

    print('RF: {}, CoCoNet: {}'.format(ml_score, coco_score))

    

if __name__ == '__main__':
    main()
