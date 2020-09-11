'''
Data loader object to aid plotting
'''

from pathlib import Path

from Bio import SeqIO
import numpy as np
import pandas as pd
import sklearn.metrics


METRICS = ['adjusted_rand_score', 'homogeneity_score', 'completeness_score']

class Loader:
    '''
    Data loader for clustering results
    Also loads Maxbin2, CONCOCT and Metabat2 results
    '''

    def __init__(self, paths=None):
        self.paths = paths

    def is_sim(self):
        '''
        Determine if it's a sim from the folder name
        '''
        return self.name.replace('_', '').isdigit()

    def get_sim_info(self):
        '''
        Retrieve sim parameters
        '''

        if self.is_sim():
            meta_names = ['bins', 'coverage', 'samples', 'iter']
            meta_vals = self.name.split('_')
            return (meta_names, meta_vals)
        return ([], [])

    def load_cluster_results(self, bins):
        '''
        Retrieves clustering results
        '''

        results = pd.read_csv(self.cfg.io['assignments'], header=None, names=['contigs', 'clusters'])

        self.data['pred'] = results.clusters
        self.data['truth'] = pd.Series({ctg: ctg.split('|')[0] for ctg in results.clusters.index})

    def load_truth(self):
        '''
        Finds true labels for:
        - simulation / true data
        - nn / clustering
        '''

        if not self.is_sim():
            truth_file = Path('{}/truth.csv'.format(self.cfg.io['output']))

            if not truth_file.is_file():
                print('Missing ground truth at {}'.format(truth_file))
                raise FileNotFoundError

            truth = pd.read_csv(truth_file, header=None, names=['contigs', 'clusters']).set_index('contigs').clusters
            common_idx = np.intersect1d(truth.index, self.data['pred'].index)
            truth = truth.loc[common_idx]

            # Check if the contig names are consistent between "truth" and the other methods results
            for method in self.data:
                if len(np.intersect1d(common_idx, self.data[method].index)) < len(common_idx):
                    print("Something went wrong: missing ids in {}".format(method))
                self.data[method] = self.data[method].loc[common_idx]
        else:
            truth = self.data['pred'].index.str.split('|').str[0]

        self.data['truth'] = pd.Series(pd.factorize(truth)[0], index=self.data['pred'].index)

    def load_concoct(self):
        '''
        Load CONCOCT binning results
        '''

        filepath = Path('{}/concurrence_results/{}/CONCOCT/concoct_clustering_gt1000.csv'
                        .format(ROOT_DIR, self.name))
        assignments = pd.read_csv(filepath)
        assignments.columns = ['contigs', 'clusters']

        self.data['CONCOCT'] = assignments.set_index('contigs').clusters

    def load_metabat2(self):
        '''
        Load Metabat2 binning results
        '''

        res_dir = Path('{}/concurrence_results/{}/Metabat2'.format(ROOT_DIR, self.name))

        assignments = []
        for i, fasta in enumerate(res_dir.glob('*.fa')):
            for seq in SeqIO.parse(fasta, 'fasta'):
                assignments.append([seq.id, i])

        if not assignments:
            assignments = np.vstack((self.data['pred'].index.values,
                                     np.arange(len(self.data['pred'])))).T

        else:
            assignments = np.array(assignments)
            cluster_imax = assignments[:, 1].astype(int).max() + 1
            ctg_rest = np.setdiff1d(self.data['pred'].index.values, assignments[:, 0])
            rest = np.array([[ctg, cluster_imax + j] for (j, ctg) in enumerate(ctg_rest)])
            assignments = np.vstack((assignments, rest))

        assignments = pd.DataFrame(assignments, columns=['contigs', 'clusters'])
        self.data['Metabat2'] = assignments.set_index('contigs').clusters

    def load_all(self, mode):
        '''
        Puts together multiple loading instructions
        '''

        if mode == 'nn':
            self.load_nn_results()
        else:
            self.load_cluster_results()
            self.load_concoct()
            self.load_metabat2()
            self.load_truth()

        ctg_order = self.data['truth'].index
        for method in self.data:
            self.data[method] = self.data[method].loc[ctg_order]

    def get_metrics(self, mode):
        '''
        Compute metrics for clustering or nn
        '''

        metrics = []
        (meta_names, meta_vals) = self.get_meta()

        for method in self.data:
            if method == 'truth':
                continue

            values = [
                [method, metric, getattr(sklearn.metrics, metric)(self.data['truth'], self.data[method])] + meta_vals
                for metric in METRICS[mode]
            ]

            metrics += values

        metrics = pd.DataFrame(metrics, columns=['method', 'metric', 'score'] + meta_names)

        return metrics

    def get_nb_clusters(self):
        
        res = {}
        for method in self.data:
            cluster_info = self.data[method].value_counts()
            res[method] = cluster_info.agg(
                {'number_of_bins': len,
                 'average_bin_size': 'mean',
                 'max_bin_size': 'max'})

        return pd.DataFrame(res)
