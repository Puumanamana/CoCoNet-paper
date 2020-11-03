import os
from pathlib import Path

from sklearn.metrics import adjusted_rand_score, homogeneity_score, completeness_score
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
    coconet.learn(cfg)
    coconet.precompute_latent_repr(cfg)

    return cfg

def compute_scores(bins_file):
    assignments = pd.read_csv(bins_file, header=None, names=['contig', 'pred'])
    assignments['truth'] = pd.factorize(assignments.contig.str.rpartition('|')[0])[0]

    return dict(
        ARI=adjusted_rand_score(assignments.truth, assignments.pred),
        homogeneity=homogeneity_score(assignments.truth, assignments.pred),
        completeness=completeness_score(assignments.truth, assignments.pred)
    )

def compute_SA_scores(bins_file, truth_file):
    assignments = pd.read_csv(bins_file, header=None, names=['contig', 'pred'])
    reference = pd.read_csv(truth_file, header=None, names=['contig', 'truth'])
    reference['truth'] = pd.factorize(reference.truth)[0]

    df = assignments.merge(reference, how='inner')

    return dict(
        ARI=adjusted_rand_score(df.truth, df.pred),
        homogeneity=homogeneity_score(df.truth, df.pred),
        completeness=completeness_score(df.truth, df.pred)
    )

def main():
    '''
    '''

    args = parser.parse_args()

    os.environ['COCONET_CONTINUE'] = 'Y'
    config = coconet_init(args)
    os.environ['COCONET_CONTINUE'] = 'N'

    sherpa_outdir = Path(args.output, 'sherpa')
    sherpa_outdir.mkdir(exist_ok=True)

    parameters = [
        sherpa.Ordinal('theta', list(np.arange(0.1, 1, 0.1))),
        sherpa.Ordinal('gamma1', list(np.arange(0.1, 1, 0.1))),
        sherpa.Ordinal('gamma2', list(np.arange(0.1, 1, 0.1))),
        sherpa.Ordinal('max_neighbors', list(50 + np.arange(0, 500, 50))),
    ]

    if args.vote_threshold is not None:
        parameters.append(sherpa.Continuous('vote_threshold', [0.01, 0.99]))
    
    algorithm = sherpa.algorithms.RandomSearch()

    study = sherpa.Study(parameters=parameters,
                         algorithm=algorithm,
                         lower_is_better=False)

    params = ['gamma1', 'gamma2', 'theta', 'max_neighbors'] + ['vote_threshold']*(args.vote_threshold is not None)
    
    for i, trial in enumerate(study):
        if trial.parameters['gamma1'] >= trial.parameters['gamma2']:
            continue

        # Update hyperparameters
        print(f"\nIteration #{i}")

        for param in params:
            setattr(config, param, trial.parameters[param])
            print(f'{param}={getattr(config, param):.1g}')

        # Run the algorithm
        coconet.cluster(config)

        # Compute the score
        if 'aloha' in str(args.fasta).lower():
            truth_file = Path(Path(args.fasta).parent, 'truth.csv')
            validation_scores = compute_SA_scores(config.io['assignments'],
                                                  str(truth_file))
        else:
            validation_scores = compute_scores(config.io['assignments'])

        print(', '.join(f'{m}={v:.3g}' for (m, v) in validation_scores.items()))

        # Add to study
        study.add_observation(trial=trial, objective=validation_scores.pop('ARI'),
                              context=validation_scores)
        study.finalize(trial)
        study.save(f'{str(sherpa_outdir)}')
    
if __name__ == '__main__':
    main()
