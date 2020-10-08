# Run the CAMISIM metagenomics simulation pipeline

## Pre-requisites
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Singularity](https://singularity.lbl.gov)

## Running the test

```bash
nextflow run nf-camisim -profile singularity \
    --outdir camisim-test
```

## Running the simulation with paper settings

```bash
nextflow run nf-camisim -profile singularity \
    --outdir output-camisim-pipeline \
    --n_genomes 500,2000 \
    --n_samples 3,15 \
    --x_coverage 4,10 \
    --n_replicates 10
```