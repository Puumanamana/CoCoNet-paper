from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt
import seaborn as sns

from coconet import coconet, parser
from coconet.core.config import Configuration
from coconet.log import setup_logger


sns.set(style='ticks', font_scale=0.7)

def coconet_init(args):
    setup_logger('<CoCoNet>', Path(args.output, 'CoCoNet.log'), args.loglvl)
    cfg = Configuration()
    cfg.init_config(**vars(args))
    cfg.to_yaml()

    coconet.preprocess(cfg)
    coconet.make_train_test(cfg)
    coconet.learn(cfg)
    coconet.precompute_latent_repr(cfg)

    return cfg

def compute_distances(feature):
    handle = feature.get_handle('latent')
    data = np.stack([np.array(handle.get(ctg)[:]) for ctg in handle.keys()])

    contig_centers = np.mean(data, axis=1)
    pairwise_distances = euclidean_distances(contig_centers)

    return pairwise_distances

def plot_distance_vs_prob(distances, truth):

    n_ctg = truth.shape[0]

    mask = np.greater(*np.meshgrid(range(n_ctg), range(n_ctg)))
    data = pd.DataFrame({name: dist[mask] for name, dist in distances.items()})
    data.columns = [f'{col.title()} space' for col in data.columns]
    data['truth'] = ['Same bin' if x else 'Different bin' for x in truth[mask]]

    g = sns.pairplot(data, hue='truth', diag_kind='hist',
                     plot_kws=dict(s=5, linewidth=0, alpha=0.8))
    g.set(xscale='log', yscale='log')
    g._legend.set_title('')
    
    g.savefig('latent_space.png', dpi=300)
    plt.show()

def main():

    args = parser.parse_args()

    cfg = coconet_init(args)

    features = cfg.get_features()

    distances = {feature.name: compute_distances(feature) for feature in features}

    contig_names = features[0].get_contigs()
    true_genomes = np.array([x.split('|')[0] for x in contig_names])
    n_ctg = len(contig_names)

    truth = np.tile(true_genomes, (n_ctg, 1)) == np.tile(true_genomes.reshape(-1, 1), (1, n_ctg))
    plot_distance_vs_prob(distances, truth)
        
if __name__ == '__main__':
    main()
