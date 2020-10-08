# Comparing binning outcome on CoCoNet paper

## Pre-requisites
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Singularity](https://singularity.lbl.gov)

## Comparing results on simulated data

```bash
nextflow run nf-binning/binning-workflow.nf -profile singularity \
     -with-report \
     -entry sim \
     --folder data/simulations
```

## Comparing results on station-aloha data

```bash
nextflow run nf-binning/binning-workflow.nf -profile singularity
     -with-report \
     -entry single \
     --folder data/station-ALOHA
```