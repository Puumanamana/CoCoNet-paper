# Repository to reproduce the figures of the CoCoNet paper

## Pre-requisites

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Singularity](https://singularity.lbl.gov)

The rest of the dependencies are already available on [DockerHub](https://hub.docker.com). The containers are defined in the "environment" folder.

For each figure, you will need to change to the corresponding directory before running the commands. For example, to reproduce the simulations, you need to do `cd camisim-simulation` first.

Finally, for all remaning scripts, you will need python >=3.6 with the following packages:
- [numpy](https://anaconda.org/conda-forge/numpy)
- [pandas](https://anaconda.org/conda-forge/pandas)
- [h5py](https://anaconda.org/conda-forge/h5py)
- [bipython](https://anaconda.org/conda-forge/biopython)
- [scikit-learn](https://anaconda.org/conda-forge/scikit-learn)
- [pingouin](https://anaconda.org/conda-forge/pingouin)
- [blast](https://anaconda.org/bioconda/blast)
- [parameters-sherpa](https://pypi.org/project/parameter-sherpa)
- [matplotlib](https://anaconda.org/conda-forge/matplotlib)
- [seaborn](https://anaconda.org/conda-forge/seaborn)
- [pytorch](https://pypi.org/project/parameter-sherpa/)
- [CoCoNet](https://pypi.org/project/coconet-binning)

## Run the simulations

Folder: *data-collection/camisim-simulation*

Simulation for binning accuracy between methods:

```bash
nextflow run nf-camisim -profile singularity \
         --outdir camisim-all \
         --n_genomes 500,2000 \
         --n_samples 4,15 \
         --x_coverage 3,10 \
         --n_replicates 10
```

Supplementary figure S8 can be directly generated with the following command:

```bash
./camisim-bin-size-distr.py --metadata camisim-all/camisim_*/metadata.csv
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

Folder: *data-collection/station-ALOHA-preprocessing*

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

For supplementary figure S1 (fragment length distribution), you can run the following command:

```bash
make tlen
```

## NCBI viral RefSeq database (for Figures 1, S7 and S9)

Folder: *data-collection/viral-refseq*

Download the database from NCBI
```bash
./download_db.sh
```

Remove any N nucleotides et replaces any other ambiguous IUPAC nucleotide at random among the available choices
```bash
./clean_db.py --fasta refseq-viral.fasta --output refseq-viral-ACGT.fasta
```

## Figure 1

Folder: *kmer-distance-distribution*

```bash
./plot_composition_separation.py --db <path to refseq db>
```

## Figures 5-7

Folder: *binning-comparison*

Running the binning on the simulation
```bash
nextflow run nf-binning -profile singularity \
         -entry sim \
         -with-report \
         --root <path to simulation folder> \
         --outdir binning-on-sim
```

## Figure 5 and Supplementary Figure S2 (FN/FP/TN/TP):
```bash
scripts/plot_nn.py --results binning-on-sim/assessment/nn_scores.csv
```

## Figure 6:
```bash
scripts/plot_clustering.py --results binning-on-sim/assessment/clustering_scores.csv --sim
```

Running the binning on the station ALOHA dataset
```bash
nextflow run nf-binning -profile singularity \
         -entry single \
         -with-report \
         --with_checkv \
         --fasta <path to station ALOHA data>/*.fasta \
         --coverage "<path to station ALOHA data>/coverage/*.bam" \
         --truth <path to station ALOHA data>/mapping_illumina-vs-ONT/mapping.csv \
         --custom_preproc "--tlen-range 200 500" \
         --outdir binning-on-SA
```

## Figure 7A and 7B:
```bash
scripts/plot_clustering.py --results binning-on-SA/assessment/clustering_scores.csv
scripts/plot_checkv binning-on-SA/assessment/checkV
```

## Supplementary Figure S3 (effect of coverage variability)

Folder: *coverage-variability-effect*

Generate the results:
```bash
./runner.sh
```

Plot:
```bash
./plot.py results/scores.csv
```

## Supplementary Figure S4 (neighbors in latent space are in the same bin)

Folder: *latent-space-neighbors*

```bash
./main.py \
    --fasta <path to one simulation folder>/assembly.fasta \
    --h5 $<path to one simulation folder>/coverage_contigs.h5
```

## Supplementary Figure S5 (effect of contig length on binning accuracy)
```bash
scripts/plot_effect_of_contig_length.py --root-dir binning_on_sim
```

## Supplementary Figure S6 (hyperparameters optimization)

Folder: *hyperparameters*
```bash
make sim INPUT=/path/to/camisim/simulation/folder
./plot_heatmaps.py --data <path to hyperparameter optimization result table> --plot-all
```

## Supplementary Figure S7 (Krona taxonomy barplot)

Folder: None

This figure requires some additional tools:
- [krona](https://anaconda.org/bioconda/krona)
- [entrez-direct](https://anaconda.org/bioconda/entrez-direct)

Extract taxids:
```bash
grep '^>' <path to refseq fasta> | cut -d ' ' -f1 | cut -c2- \
    | epost -db nuccore \
    | esummary -db nuccore \
    | xtract -pattern DocumentSummary -element Caption,TaxId \
    > taxids.tsv 
```

Plot figure with Krona:
```bash
ktImportTaxonomy taxids.tsv  
```

## Supplementary Figure S9 (speed/memory)

Folder: *speed-memory*

Generate the simulations with 7 threads:
```bash
make -j 7
```

Process the simulations
```bash
nextflow run ../binning-comparison/nf-binning -profile singularity \
    -entry sim \
    -with-report \
    --root simulations
```

Nextflow will generate an html report, which includes a summary table. This will need to be converted to csv format.
Plot:
```bash
./plot-speed-mem.py <path to result table>
```
