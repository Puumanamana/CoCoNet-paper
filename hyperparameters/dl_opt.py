import os
from math import ceil
from pathlib import Path

from sklearn.metrics import accuracy_score, confusion_matrix
import pandas as pd
import numpy as np

import sherpa

from coconet import coconet, parser
from coconet.core.config import Configuration
from coconet.log import setup_logger


def coconet_init(args):
    setup_logger('CoCoNet', Path(args.output, 'CoCoNet.log'), args.loglvl)
    cfg = Configuration()
    cfg.init_config(**vars(args))
    cfg.to_yaml()

    coconet.preprocess(cfg)
    coconet.make_train_test(cfg)

    return cfg

def compute_scores(pred_file):
    res = pd.read_csv(pred_file)
    pred_bin = (res.combined > 0.5).astype(int)
    matrix = confusion_matrix(res.truth, pred_bin)
    
    return dict(
        acc=accuracy_score(res.truth, pred_bin),
        TP=matrix[1, 1],
        TN=matrix[0, 0],
        FN=matrix[1, 0],
        FP=matrix[0, 1],                 
    )


def main():
    '''
    '''

    args = parser.parse_args()
    
    os.environ['COCONET_CONTINUE'] = 'Y'
    config = coconet_init(args)
    os.environ['COCONET_CONTINUE'] = 'N'

    sherpa_outdir = Path(args.output, 'dl-sherpa')
    sherpa_outdir.mkdir(exist_ok=True)

    parameters = [
        sherpa.Ordinal('compo1-factor', [2, 4, 8]),
        sherpa.Ordinal('compo2', [8, 16, 32, 64]),
        sherpa.Ordinal('cover-kernel', [2, 4, 8]),
        sherpa.Ordinal('cover-stride-ratio', [0.5, 1]),
        sherpa.Ordinal('cover-filters', [16, 32, 64]),
        # sherpa.Ordinal('wsize', [2, 4, 8, 16]),
        # sherpa.Continuous('wstep', [0.2, 1], scale='log'),                        
        sherpa.Ordinal('merge', [8, 16, 32]),
        sherpa.Discrete('kmer', [4, 5]),
        # sherpa.Continuous('lr', [1e-4, 1e-3], scale='log'),
        sherpa.Ordinal('bs', [128, 256, 512]),
    ]
    
    algorithm = sherpa.algorithms.GPyOpt()

    study = sherpa.Study(parameters=parameters,
                         algorithm=algorithm,
                         lower_is_better=False)
    
    for i, trial in enumerate(study):
        # Update hyperparameters
        print(f"\nIteration #{i}")
        
        config.compo_neurons = np.array([int(trial.parameters['compo1-factor']), 1]) * int(trial.parameters['compo2'])
        config.cover_neurons = config.compo_neurons
        config.merge_neurons = int(trial.parameters['merge'])
        config.cover_filters = int(trial.parameters['cover-filters'])
        config.cover_kernel = int(trial.parameters['cover-kernel'])
        config.cover_stride = ceil(trial.parameters['cover-stride-ratio'] * config.cover_kernel)
        # config.wsize = trial.parameters['wsize']
        # config.wstep = ceil(trial.parameters['wstep'] * config.wsize)
        config.kmer = int(trial.parameters['kmer'])
        config.batch_size = int(trial.parameters['bs'])
        config.learning_rate = 1e-3 # trial.parameters['lr']        

        # Run the algorithm
        coconet.learn(config)

        validation_scores = compute_scores(config.io['nn_test'])

        print(', '.join(f'{m}={v:.2%}' if m == 'acc' else f'{m}={v:,}'
                        for (m, v) in validation_scores.items()))
        
        # Add to study
        study.add_observation(trial=trial, objective=validation_scores.pop('acc'),
                              context=validation_scores)
        study.finalize(trial)
        study.save(f'{str(sherpa_outdir)}')

        
if __name__ == '__main__':
    main()
