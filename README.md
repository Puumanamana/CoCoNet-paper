# Repository to reproduce the figures of the CoCoNet paper

## Pre-requisites

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Singularity](https://singularity.lbl.gov)

The rest of the dependencies are already available on [DockerHub](https://hub.docker.com). The containers are defined in the "environment" folder.

For each figure, you will need to change to the corresponding directory before running the commands. For example, to reproduce the simulations, you need to do `cd camisim-simulation` first.

## Run the simulations

Folder: *camisim-simulation*

Simulation for binning accuracy between methods:

```bash
nextflow run nf-camisim -profile singularity \
         --outdir camisim-all \
         --n_genomes 500,2000 \
         --n_samples 4,15 \
         --x_coverage 3,10 \
         --n_replicates 10
```

Simulation for hyperparameter optimization:

```bash
nextflow run nf-camisim -profile singularity \
         --outdir camisim-for-opt \
         --n_genomes 4000 \
         --n_samples 5 \
         --x_coverage 6 \
         --n_replicates 1
```

## Download and preprocess Station ALOHA dataset

Folder: *station-ALOHA-preprocessing*

The nf-station-aloha pipeline will:
    - Download both the illumina
    - Trim and assemble the illumina reads
    - Align the trimmed reads on the assembly
    - Download and de-duplicate the ONT assembly
    - Align the illumina contigs against the ONT assembly
    - Create a mapping file ("true" contig assignments)
    
```bash
nextflow run nf-station-aloha --run-assembly --outdir station-ALOHA -profile singularity
```

For the fragment length distribution in the supplementary, you can run the following command:

```bash
make tlen
```

## Figure 1

Folder: *kmer-distance-distribution*

Download the database. If any backlashes are automatically added before the curly braces, you will need to remove them.

```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{1,2}.1.genomic.fna.gz
cat viral*.gz > refseq-viral.fasta && rm -f viral*.gz
```

Clean the database by removing N nucleotides and converting any other ambiguous nucleotides 
```bash
python clean_db.py --fasta refseq-viral.fasta --output refseq-viral-ACGT.fasta
```

Run the main script
```bash
python plot_composition_separation.py --db refseq-viral-ACGT.fasta
```

## Figures 5-7

Folder: *binning-comparison*

First, you will need to create a "data" folder and symlink the output folders of both the simulation and station ALOHA in it.


Running the binning on the simulation
```bash
nextflow run nf-binning -profile singularity \
         -entry sim \
         -with-report \
         --root data/simulations \
         --outdir binning-on-sim
```

Figure 5 and supplementary figure with FN/FP/TN/TP:
```bash
python scripts/plot_nn.py --results binning-on-sim/assessment/nn_scores.csv
```

Figure 6:
```bash
python scripts/plot_clustering.py --results binning-on-sim/assessment/clustering_scores.csv --sim
```

Running the binning on the station ALOHA dataset
```bash
nextflow run nf-binning -profile singularity \
         -entry single \
         -with-report \
         --with_checkv \
         --fasta data/station-ALOHA/*.fasta \
         --coverage "data/station-ALOHA/*.bam" \
         --custom_preproc "--tlen-range 200 500" \
         --outdir binning-on-SA
```

Figure 7A and 7B:
```bash
python scripts/plot_clustering.py --results binning-on-SA/assessment/clustering_scores.csv
python scripts/plot_checkv binning-on-SA/assessment/checkV
```

## Supplementary (neighbors in latent space are in the same bin)

Folder: *latent-space-neighbors*

```bash
make INPUT=/path/to/one/simulation/folder
```

## Supplementary (effect of coverage variability)

Folder: *coverage-variability-effect*

```bash
bash runner.sh /path/to/simulations/root/folder
python plot.py
```

## Supplementary (effect of contig length on binning accuracy)

Folder: *?*

```bash

```

