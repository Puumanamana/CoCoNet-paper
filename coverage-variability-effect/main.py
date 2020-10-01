from pathlib import Path
import logging

import numpy as np
import pandas as pd
import torch
import torch.optim as optim
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    confusion_matrix,
    roc_auc_score
)
from coconet import coconet, parser
from coconet.core.config import Configuration
from coconet.core.generators import CoverageGenerator
from coconet import dl_util
from coconet.log import setup_logger

TEST_FREQ = 100


def coconet_init(args):
    setup_logger('CoCoNet', Path(args.output, 'CoCoNet.log'), args.loglvl)
    cfg = Configuration()
    cfg.init_config(**vars(args))
    cfg.to_yaml()

    coconet.preprocess(cfg)
    coconet.make_train_test(cfg)

    return cfg

def main():
    '''
    '''

    args = parser.parse_args()

    cfg = coconet_init(args)
    logger = logging.getLogger('learning')
    
    (x_test, x_train) = (get_generator(cfg, mode) for mode in ['test', 'train'])
    (y_test, y_train) = (get_truth(cfg.io['pairs'][mode]) for mode in ['test', 'train'])

    x_test = dict(with_var=next(x_test))
    x_test['no_var'] = [torch.ones_like(xi) * xi.mean(axis=2).unsqueeze(-1)
                        for xi in x_test['with_var']]
    y_test = y_test.numpy().astype(int)
    
    scores = []
    models = dict()
    optimizers = dict()

    n_batch = 1 + dl_util.get_npy_lines(cfg.io['pairs']['train']) // cfg.batch_size

    for coverage_mode in x_test:
        models[coverage_mode] = dl_util.initialize_model(
            'coverage', cfg.get_input_shapes(), cfg.get_architecture()
        ).train()
        optimizers[coverage_mode] = optim.Adam(models[coverage_mode].parameters(),
                                               lr=cfg.learning_rate)

    for i, batch_x in enumerate(x_train, 1):
        for coverage_mode in x_test:
            optimizers[coverage_mode].zero_grad()

            truth_b = y_train[(i-1)*cfg.batch_size:i*cfg.batch_size]

            if coverage_mode == 'with_var':
                pred_b = models[coverage_mode](*batch_x)
            else:
                template = torch.ones_like(batch_x[0])
                batch_x_mean = [template*xi.mean(axis=2).unsqueeze(-1) for xi in batch_x]
                pred_b = models[coverage_mode](*batch_x_mean)

            loss = models[coverage_mode].compute_loss(pred_b, truth_b)
            loss.backward()
            
            optimizers[coverage_mode].step()

        # Get test results
        if (i % TEST_FREQ == 0 and i > 0):
            for coverage_mode in models:
                pred = dl_util.run_test(models[coverage_mode], x_test[coverage_mode])
                metrics = get_scores(y_test, pred['coverage'])
                metrics.update(mode=coverage_mode, batch=i)
                scores.append(metrics)

                logger.info((
                    f'(Batch #{i:,} / {n_batch:,}, {coverage_mode}) '
                    f'accuracy={metrics["accuracy"]:.2%}, AUC={metrics["AUC"]:.2%}'
                ))

    scores = pd.DataFrame(scores)
    scores.to_csv(f'results/test-{cfg.io["fasta"].parent.stem}.csv', index=False)

def get_truth(pairs_file):
    ctg_names = np.load(pairs_file)['sp']

    sim = all('|' in ctg for ctg in ctg_names[:, 0])

    if sim:
        truth = np.zeros(len(ctg_names))
        for i, (ctg1, ctg2) in enumerate(ctg_names):
            truth[i] = ctg1 == ctg2
    else:
        truth = ctg_names[:, 0] == ctg_names[:, 1]

    return torch.from_numpy(truth.astype(np.float32)[:, None])

def get_generator(cfg, mode):
    pairs = cfg.io['pairs'][mode]

    batch_size = cfg.batch_size
    if mode == 'test':
        batch_size = np.load(pairs).shape[0]
        
    return CoverageGenerator(
        pairs,
        cfg.io['h5'],
        batch_size=batch_size,
        load_batch=cfg.load_batch,
        wsize=cfg.wsize, wstep=cfg.wstep
    )

def get_scores(truth, pred):
    pred_bin = pred.round()
    summary = confusion_matrix(truth, pred_bin)

    return dict(
        accuracy=accuracy_score(truth, pred_bin),
        AUC=roc_auc_score(truth, pred),
        precision=precision_score(truth, pred_bin),
        recall=recall_score(truth, pred_bin),
        TN=summary[0, 0],
        TP=summary[1, 1],
        FN=summary[1, 0],
        FP=summary[0, 1],
    )


if __name__ == '__main__':
    main()
