from pathlib import Path
from sklearn.metrics import adjusted_rand_score
import pandas as pd
import sherpa

from coconet import coconet, parser
from coconet.core import config

def coconet_init(args):
    cfg = config.Configuration()
    cfg.init_config(mkdir=True, **vars(args))
    cfg.to_yaml()

    coconet.preprocess(cfg)
    coconet.make_train_test(cfg)
    coconet.learn(cfg)

    return cfg

def compute_score(bins_file):
    assignments = pd.read_csv(bins_file, header=None, names=['contig', 'pred'])
    assignments['truth'] = pd.factorize(assignments.contig.str.rpartition('|')[0])[0]
    
    return adjusted_rand_score(assignments.truth, assignments.pred)

def compute_SA_score(bins_file, truth_file):
    assignments = pd.read_csv(bins_file, header=None, names=['contig', 'pred'])
    truth = pd.read_csv(truth_file, header=None, names=['contig', 'truth'])

    df = assignments.merge(truth, how='inner')

    return adjusted_rand_score(df.truth, df.pred)

def main():
    '''
    '''

    args = parser.parse_args()

    name = Path(args.fasta).parent.stem
    config = coconet_init(args)

    parameters = [
        sherpa.Continuous('theta', [0.01, 0.99]),
        sherpa.Continuous('gamma1', [0.05, 1]),        
        sherpa.Continuous('gamma2', [0.05, 1]), # gamma2 = gamma1 + delta_g
        sherpa.Discrete('max_neighbors', [1, 500])
    ]
    algorithm = sherpa.algorithms.GPyOpt()

    study = sherpa.Study(parameters=parameters,
                         algorithm=algorithm,
                         lower_is_better=False)

    for i, trial in enumerate(study):
        # Update hyperparameters
        print(f"\nIteration #{i}")
        config.gamma1 = trial.parameters['gamma1']
        config.gamma2 = trial.parameters['gamma2']
        config.theta = trial.parameters['theta']
        config.max_neighbors = trial.parameters['max_neighbors']        
        print(f'''
        gamma1={config.gamma1:.4f}, 
        gamma2={config.gamma2:.4f}, 
        theta={config.theta:.4f}, 
        max_neighbors={config.max_neighbors}
        ''')

        # Run the algorithm
        coconet.cluster(config, force=True)

        # Compute the score
        if 'station_aloha' in str(args.fasta).lower():
            truth_file = Path(Path(args.fasta).parent, 'truth.csv')
            validation_error = compute_SA_score(config.io['assignments'],
                                                str(truth_file))
        else:
            validation_error = compute_score(config.io['assignments'])

        # Add to study
        study.add_observation(trial=trial, objective=validation_error)
        study.finalize(trial)
        study.save(f'sherpa-results-{name}')
        
if __name__ == '__main__':
    main()
