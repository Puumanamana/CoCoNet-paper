import configparser
import pandas as pd

project_path = "."

def generate_config(size=1,n_samples=5,n_genomes=6000,ds_name=None, cami_src=None):
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
        'readsim': '{}/tools/art_illumina-2.3.6/art_illumina'.format(cami_src),
        'error_profiles': '{}/tools/art_illumina-2.3.6/profiles'.format(cami_src),
        'samtools': '{}/tools/samtools-1.3/samtools'.format(cami_src),
        'profile': 'mbarc',
        'size': size,
        'type': 'art',
        'fragments_size_mean': 400,
        'fragment_size_standard_deviation': 10
    }
    config['CommunityDesign'] = {
        'ncbi_taxdump': '{}/tools/ncbi-taxonomy_20190708.tar.gz'.format(cami_src),
        'number_of_samples': n_samples
    }
    config['community0'] = {
        'metadata': '{}/metadata_camisim.tsv'.format(project_path),
        'id_to_genome_file': '{}/id_to_genome.tsv'.format(project_path),
        'genomes_total': n_genomes,
        'genomes_real': n_genomes,
        'max_strains_per_otu': 1,
        'ratio': 1,
        'mode': 'differential',
        'log_mu': 1,
        'log_sigma': 3
    }
    with open('config.ini', 'w') as configfile:
        config.write(configfile)

def generate_id_to_genome(genomes,out_dir="."):
    mapping = pd.Series({ f.split(".")[0]: "contigs_fasta/{}.fasta".format(f)
                          for f in genomes })
    mapping.to_csv("{}/id_to_genome.tsv".format(out_dir),sep="\t",header=False)

def generate_metadata(out_dir="."):

    meta_camisim = pd.read_csv("{}/id_to_genome.tsv".format(out_dir),
                               sep="\t",index_col=0,header=None)
    metadata = pd.DataFrame(columns=["OTU","NCBI_ID","novelty_category"],
                            index=meta_camisim.index)
    metadata.index.name = "genome_ID"

    metadata["OTU"] = range(metadata.shape[0])
    metadata["NCBI_ID"] = 10239
    metadata["novelty_category"] = 'new_species'

    metadata.to_csv("{}/metadata_camisim.tsv".format(out_dir),sep="\t")
        
if __name__ == '__main__':
    import argparse
    from Bio import SeqIO

    parser = argparse.ArgumentParser(description='Generate metadata for camisim')
    parser.add_argument('--cov_lvl', type=float)
    parser.add_argument('--n_samples', type=int)
    parser.add_argument('--n_genomes', type=int)
    parser.add_argument('--db', type=str)
    parser.add_argument('--cami_src', type=str)        
    args = parser.parse_args()

    all_genomes = pd.Series({ seq.id: len(seq.seq) for seq in SeqIO.parse(args.db,'fasta') })
    sub_genomes = all_genomes.sample(args.n_genomes)
                            
    ds_name = 'sim_ng={:g}_cov={:g}X_nsamples={:g}'.format(args.n_genomes, args.cov_lvl, args.n_samples)
    
    generate_config(size=args.cov_lvl*sub_genomes.sum()/1e9,
                    n_samples=args.n_samples,
                    n_genomes=len(sub_genomes),
                    ds_name=ds_name,
                    cami_src=args.cami_src)
    generate_id_to_genome(sub_genomes.index)
    generate_metadata()
