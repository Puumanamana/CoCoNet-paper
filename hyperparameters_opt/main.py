import os
from pathlib import Path

from sklearn.metrics import adjusted_rand_score
import pandas as pd

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

def compute_score(bins_file):
    assignments = pd.read_csv(bins_file, header=None, names=['contig', 'pred'])
    assignments['truth'] = pd.factorize(assignments.contig.str.rpartition('|')[0])[0]

    return adjusted_rand_score(assignments.truth, assignments.pred)

def compute_SA_score(bins_file, truth_file):
    assignments = pd.read_csv(bins_file, header=None, names=['contig', 'pred'])
    reference = pd.read_csv(truth_file, header=None, names=['contig', 'truth'])
    reference['truth'] = pd.factorize(reference.truth)[0]

    df = assignments.merge(reference, how='inner')

    return adjusted_rand_score(df.truth, df.pred)

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
        sherpa.Continuous('theta', [0.01, 0.99]),
        sherpa.Continuous('gamma1', [0.05, 1]),        
        sherpa.Continuous('gamma2', [0.05, 1]), # gamma2 = gamma1 + delta_g
        sherpa.Discrete('max_neighbors', [1, 500])
    ]
    if args.vote_threshold is not None:
        parameters.append(sherpa.Continuous('vote_threshold', [0.01, 0.99]))
    
    algorithm = sherpa.algorithms.GPyOpt()

    study = sherpa.Study(parameters=parameters,
                         algorithm=algorithm,
                         lower_is_better=False)

    for i, trial in enumerate(study):
        # Update hyperparameters
        print(f"\nIteration #{i}")

        for param in (['gamma1', 'gamma2', 'theta', 'max_neighbors'] +
                      ['vote_threshold']*(args.vote_threshold is not None)):
            setattr(config, param, trial.parameters[param])
            print(f'{param}={getattr(config, param):.1g}')

        # Run the algorithm
        coconet.cluster(config)

        # Compute the score
        if 'aloha' in str(args.fasta).lower():
            truth_file = Path(Path(args.fasta).parent, 'truth.csv')
            validation_score = compute_SA_score(config.io['assignments'],
                                                str(truth_file))
        else:
            validation_score = compute_score(config.io['assignments'])

        print(f'Validation score: {validation_score:.2f}')
        # Add to study
        study.add_observation(trial=trial, objective=validation_score)
        study.finalize(trial)
        study.save(f'{str(sherpa_outdir)}')

        
if __name__ == '__main__':
    main()
