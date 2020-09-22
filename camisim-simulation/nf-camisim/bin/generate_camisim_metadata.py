#!/usr/bin/env python

import configparser
import pandas as pd
import argparse
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='Generate metadata for camisim')
    parser.add_argument('--name', type=str, default='camisim')    
    parser.add_argument('--cov_lvl', type=float)
    parser.add_argument('--n_samples', type=int)
    parser.add_argument('--n_genomes', type=int)
    parser.add_argument('--db', type=str)
    parser.add_argument('--cami-data', type=str, default='.')
    parser.add_argument('--log-mu', type=int, default=1)
    parser.add_argument('--log-sigma', type=int, default=2)
    
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    all_genomes = pd.Series({ seq.id: len(seq.seq) for seq in SeqIO.parse(args.db,'fasta') })
    sub_genomes = all_genomes.sample(args.n_genomes)

    generate_config(size=args.cov_lvl*sub_genomes.sum()/1e9,
                    n_samples=args.n_samples,
                    n_genomes=len(sub_genomes),
                    ds_name=args.name,
                    cami_data=args.cami_data,
                    log_mu=args.log_mu,
                    log_signma=args.log_sigma)

    generate_id_to_genome(sub_genomes.index)
    generate_metadata()

def generate_config(ds_name=None, cami_data='.', project_path='.',
                    size=1, n_samples=5, n_genomes=6000,
                    log_mu=1, log_sigma=2):
    config = configparser.ConfigParser()
    config['Main'] = {
        'seed': 42,
        'phase': 0,
        'max_processors': 50,
        'dataset_id': ds_name,
        'output_directory': ds_name,
        'temp_directory': '/tmp/',
        'gsa': False,
        'pooled_gsa': True,
        'anonymous': False,
        'compress': 0
    }
    config['ReadSimulator'] = {
        'readsim': f'{cami_data}/tools/art_illumina-2.3.6/art_illumina',
        'error_profiles': f'{cami_data}/tools/art_illumina-2.3.6/profiles',
        'samtools': f'{cami_data}/tools/samtools-1.3/samtools',
        'profile': 'mbarc',
        'size': size,
        'type': 'art',
        'fragments_size_mean': 400,
        'fragment_size_standard_deviation': 10
    }
    config['CommunityDesign'] = {
        'ncbi_taxdump': f'{cami_data}/tools/ncbi-taxonomy_20190708.tar.gz',
        'number_of_samples': n_samples
    }
    config['community0'] = {
        'metadata': f'{project_path}/metadata_camisim.tsv',
        'id_to_genome_file': f'{project_path}/id_to_genome.tsv',
        'genomes_total': n_genomes,
        'genomes_real': n_genomes,
        'max_strains_per_otu': 1,
        'ratio': 1,
        'mode': 'differential',
        'log_mu': log_mu,
        'log_sigma': log_sigma
    }
    with open('config.ini', 'w') as configfile:
        config.write(configfile)

def generate_id_to_genome(genomes,out_dir="."):
    mapping = pd.Series({f.split(".")[0]: f"contigs_fasta/{f}.fasta" for f in genomes})
    mapping.to_csv(f"{out_dir}/id_to_genome.tsv", sep="\t", header=False)

def generate_metadata(out_dir="."):

    meta_camisim = pd.read_csv(f"{out_dir}/id_to_genome.tsv", sep="\t", index_col=0, header=None)
    metadata = pd.DataFrame(columns=["OTU", "NCBI_ID", "novelty_category"], index=meta_camisim.index)
    metadata.index.name = "genome_ID"

    metadata["OTU"] = range(metadata.shape[0])
    metadata["NCBI_ID"] = 10239
    metadata["novelty_category"] = 'new_species'

    metadata.to_csv(f"{out_dir}/metadata_camisim.tsv", sep="\t")


if __name__ == '__main__':
    main()
