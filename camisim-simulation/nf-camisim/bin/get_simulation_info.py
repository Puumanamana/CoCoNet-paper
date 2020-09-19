import os,sys

import h5py
import numpy as np
from progressbar import progressbar

PARENT_DIR = os.path.join(sys.path[0], '..')

def get_coverage(coverage, n_samples):
    '''
    Calculates the coverage per contig in the simulation
    '''

    # folder = "{}/input_data/{}_{}".format(PARENT_DIR, coverage, nsamples)
    folder = 'camisim_results/{}_{}'.format(coverage, nsamples)

    h5file = h5py.File('{}/coverage_contigs.h5'.format(folder), 'r')
    contigs = list(h5file.keys())


    try:
        n_samples_sim = h5file.get(contigs[0]).shape[0]

        total_len = np.zeros(n_samples_sim)
        total_bp = np.zeros(n_samples_sim)
        total_coverage = np.zeros((len(h5file), n_samples_sim))
        
        for i, ctg in progressbar(enumerate(h5file), max_value=len(h5file)):
            data = h5file.get(ctg)
            ctg_len = data.shape[1]
            ctg_cov = data[:].sum(axis=1)

            total_coverage[i] = ctg_cov / ctg_len
            total_len += ctg_len
            total_bp += ctg_cov
            
    except IndexError:
        total_coverage = np.zeros((2,2))
            
    print('''
    In theory: n_samples = {} - coverage = {}
    Observed: n_samples = {}
    Depth = {}X
    coverage = {}'''
          .format(n_samples, coverage,
                  total_coverage.shape[1], 'X '.join((total_bp/total_len).astype(str))+'X',
                  total_coverage.mean(axis=0)))
    return total_coverage

if __name__ == '__main__':
    # import argparse

    # parser = argparse.ArgumentParser()
    # parser.add_argument('--coverage', type=int)
    # parser.add_argument('--nsamples', type=int)

    # args = parser.parse_args()

    for nsamples in [3,5,10]:
        for coverage in [1,2,5,10]:
            get_coverage(coverage, nsamples)
