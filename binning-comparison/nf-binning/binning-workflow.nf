#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { REFORMAT_COVERAGE; BINS_TO_FASTA; COMPUTE_SCORES } from './util/process'
include { COCONET_PREPROCESS } from './modules/coconet/preprocess/process'
include { COCONET_RUN } from './modules/coconet/run/process' addParams(outdir: "$params.outdir/binning")
include { CONCOCT } from './modules/concoct/process' addParams(outdir: "$params.outdir/binning")
include { METABAT2 } from './modules/metabat2/process' addParams(outdir: "$params.outdir/binning")
include { MAXBIN2 } from './modules/maxbin2/process' addParams(outdir: "$params.outdir/binning")


workflow binning {
    take:
    folder

    main:
    ds = "${folder.getName()}"
    fasta = Channel.fromPath("${folder}/*.fasta").map{[ [id: ds], it]}
    bams = Channel.fromPath("${folder}/*.bam").collect().map{[ [id: ds], it]}
    h5 = Channel.fromPath("${folder}/*.{h5,hdf5}").map{[ [id: ds], it]}
    truth = file("${folder}/truth.csv")

    preproc_args = ["--min-ctg-len", "$params.min_ctg_len",
                    "--min-mapping-quality", "$params.min_mapq",
                    "--min-aln-coverage", "$params.min_map_id",
                    "--tlen-range", "$params.fl_range"]

    coconet_args = preproc_args + ["--min-prevalence", "${params.coconet.min_prevalence}",
                                   "--theta", "${params.coconet.theta}",
                                   "--gamma1", "${params.coconet.gamma1}",
                                   "--gamma2", "${params.coconet.gamma2}",
                                   "--features", "${params.coconet.features}"]
   
    // inputs preprocessing
    inputs = COCONET_PREPROCESS(
        fasta.join(bams.mix(h5)),
        [publish_dir: 'preprocessing',
         args: "${preproc_args.join(' ')} --min-prevalence 0"]
    )
    coverage = REFORMAT_COVERAGE(inputs.fasta.join(inputs.h5))

    // binning
    coconet_bins = COCONET_RUN(
        inputs.fasta.join(inputs.h5),
        [publish_dir: 'CoCoNet',
         args: coconet_args.join(' ')]
    ).bins
    
    metabat2_bins = METABAT2(
        inputs.fasta.join(coverage.metabat2),
        [publish_dir: 'Metabat2',
         args: "--minContig ${params.min_ctg_len}"]
    ).bins
    concoct_bins = CONCOCT(
        inputs.fasta.join(coverage.concoct),
        [publish_dir: 'CONCOCT',
         args: "--clusters ${params.concoct.max_clusters} --length_threshold ${params.min_ctg_len}"]
    ).bins
    maxbin2_bins = MAXBIN2(
        inputs.fasta.join(coverage.maxbin2),
        [publish_dir: 'MaxBin2',
         args: "-min_contig_length ${params.min_ctg_len}"]
    ).bins

    all_bins = coconet_bins.mix(metabat2_bins, concoct_bins, maxbin2_bins)

    // compute scores
    COMPUTE_SCORES(
        ds,
        all_bins.collect{it[1]},
        truth
    )

    // reformating
    BINS_TO_FASTA(
        coconet_bins.mix(all_bins).combine(fasta, by: 0)
    )
}


workflow all_sims {
    Channel.fromPath("${params.root_dir}/*", type: 'dir') | binning
}

workflow single {
    binning(file(params.folder, checkIfExists: true))
}

workflow {
    binning(file("${baseDir}/test"))
}
