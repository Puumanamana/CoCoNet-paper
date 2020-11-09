from pathlib import Path
import re
import logging
from collections import deque

import numpy as np
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
from coconet.core.generators import CoverageGenerator, CompositionGenerator
from coconet import dl_util

TEST_FREQ = 100
METRICS = ['accuracy', 'AUC', 'TP', 'TN', 'FP', 'FN']


def coconet_init(args):
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

    # load data
    (x_test_compo, x_train_compo) = (get_composition_generator(cfg, s)
                                     for s in ['test', 'train'])
    (x_test_cover, x_train_cover) = (get_coverage_generator(cfg, s)
                                     for s in ['test', 'train'])
    (y_test, y_train) = (get_truth(cfg.io['pairs'][s]) for s in ['test', 'train'])

    if len(cfg.features) > 1:
        x_train = zip(x_train_compo, x_train_cover)
        x_test = dict(with_var=next(zip(x_test_compo, x_test_cover)))
    else:
        x_train = x_train_cover
        x_test = dict(with_var=next(x_test_cover))

    x_test['no_var'] = flatten_coverage(x_test['with_var'])
    y_test_npy = y_test.numpy().astype(int)

    # Initialize models
    n_batch = 1 + dl_util.get_npy_lines(cfg.io['pairs']['train']) // cfg.batch_size
    nets = dict(no_var=dict(), with_var=dict())

    for mode in ['no_var', 'with_var']:
        nets[mode]['model'] = dl_util.initialize_model(
            '-'.join(cfg.features), cfg.get_input_shapes(), cfg.get_architecture()
        ).train()
        nets[mode]['optim'] = optim.Adam(nets[mode]['model'].parameters(), lr=cfg.learning_rate)
        nets[mode]['loss'] = deque(maxlen=args.patience)
        nets[mode]['over'] = False

    logger = logging.getLogger('<learning>')
    logger.info('Training started')

    for i, batch_x in enumerate(x_train, 1):

        if all(nets[mode]['over'] for mode in nets):
            break

        for mode in nets:
            if nets[mode]['over']:
                continue

            nets[mode]['optim'].zero_grad()

            truth_b = y_train[(i-1)*cfg.batch_size:i*cfg.batch_size]

            if mode == 'with_var':
                pred_b = nets['with_var']['model'](*batch_x)
            else:
                batch_x_flat = flatten_coverage(batch_x)
                pred_b = nets['no_var']['model'](*batch_x_flat)

            loss = nets[mode]['model'].compute_loss(pred_b, truth_b)
            loss.backward()

            nets[mode]['optim'].step()

        if not (
                (i % TEST_FREQ == 0) or (i == n_batch)
        ):
            continue

        # Get test results
        for mode in nets:
            if nets[mode]['over']:
                continue

            (pred, test_loss) = make_prediction(
                nets[mode]['model'], x_test[mode], y_test, args.features
            )

            losses = nets[mode]['loss']
            losses.append(test_loss)

            scores = get_scores(y_test_npy, pred)
            scores.update(mode=mode, batch=i)

            logger.info((
                f'(Batch #{i:,} / {n_batch:,}, {mode}) '
                f'accuracy={scores["accuracy"]:.2%}, AUC={scores["AUC"]:.2%}'
            ))

            if test_loss <= np.min(losses):
                nets[mode]['best_so_far'] = [i] + [scores[k] for k in METRICS]

            if (len(losses) == args.patience
                and losses[0] <= np.min(losses)):
                nets[mode]['over'] = True
                logger.info(f'{mode}: early stopping (best={nets[mode]["best_so_far"]}')

    save_best_scores(args.fasta, [(mode, nets[mode]['best_so_far']) for mode in nets],
                     folder=args.output.parent)

def save_best_scores(fasta, results, folder):
    if 'aloha' in str(fasta).lower():
        print(results)
        return

    output = Path(f'{folder}/scores.csv')
    if not output.is_file():
        header = ['genomes', 'samples', 'coverage', 'replicate',
                  'mode', 'last_batch'] + METRICS
        with open(output, 'w') as handle:
            handle.write(','.join(header) + '\n')

    with open(output, 'a') as handle:
        info = re.split(r'[_\-]', Path(fasta).parent.name)

        for (mode, scores) in results:
            data = [info[i] for i in [1, 3, 5, 6]] + [mode] + scores
            handle.write(','.join(map(str, data)) + '\n')

def flatten_coverage(inputs):
    if len(inputs[0]) == 2:
        # combined model
        (compo, cover) = inputs
    else:
        cover = inputs

    template = torch.ones_like(cover[0])
    flat_cov = [template*xi.mean(axis=2).unsqueeze(-1) for xi in cover]

    if len(inputs[0]) == 2:
        return [compo, flat_cov]

    return flat_cov
            
def make_prediction(model, x, y, ft='coverage'):

    if len(ft) == 1:
        ft = ft[0]
    else:
        ft = 'combined'

    model.eval()
    pred = model(*x)[ft]
    loss = model.loss_op(pred, y).mean().item()
    model.train()

    return (pred.detach().numpy(), loss)
    
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

def get_composition_generator(cfg, mode):
    pairs = cfg.io['pairs'][mode]

    return CompositionGenerator(
        pairs, cfg.io['filt_fasta'],
        batch_size=cfg.batch_size * (mode=='train'),
        kmer=cfg.kmer, rc=True, norm=False
    )

def get_coverage_generator(cfg, mode):
    pairs = cfg.io['pairs'][mode]

    return CoverageGenerator(
        pairs,
        cfg.io['h5'],
        batch_size=cfg.batch_size * (mode=='train'),
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
