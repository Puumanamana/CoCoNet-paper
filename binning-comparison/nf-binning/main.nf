#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { REFORMAT_COVERAGE; MERGE_BINS; COMPUTE_CLUSTERING_METRICS; COMPUTE_NN_METRICS } from './util/process'
include { COCONET_PREPROCESS } from './modules/coconet/preprocess/process'
include { COCONET_RUN } from './modules/coconet/run/process' addParams(outdir: "$params.outdir/binning")
include { CONCOCT } from './modules/concoct/process' addParams(outdir: "$params.outdir/binning")
include { METABAT2 } from './modules/metabat2/process' addParams(outdir: "$params.outdir/binning")
include { CHECKV_END_TO_END } from './modules/checkv/end_to_end/process' addParams(outdir: "$params.outdir/assessment")


workflow binning {
    take:
    fasta
    coverage
    truth

    main:
    preproc_args = ["--min-ctg-len", "$params.min_ctg_len",
                    "--min-mapping-quality", "$params.min_mapq",
                    "--min-aln-coverage", "$params.min_map_id",
                    "${params.custom_preproc}"]

    // inputs preprocessing
    inputs = COCONET_PREPROCESS(
        fasta.join(coverage),
        [publish_dir: 'preprocessing',
         args: "${preproc_args.join(' ')} --min-prevalence 0"]
    )
    coverage = REFORMAT_COVERAGE(inputs.fasta.join(inputs.h5))

    // binning
    coconet_bins = COCONET_RUN(
        inputs.fasta.join(inputs.h5),
        [publish_dir: 'CoCoNet',
         args: preproc_args.join(' ') + " ${params.coconet.args}"]
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
    
    all_bins = coconet_bins.mix(metabat2_bins, concoct_bins).groupTuple()

    // reformating
    bins_fmt = MERGE_BINS(
        all_bins.map{[[id: it[0].id], it]}.combine(inputs.fasta, by: 0).map{[it[1][0], it[1][1], it[2]]}
    )

    // compute NN metrics
    COMPUTE_NN_METRICS(
        COCONET_RUN.out.nn_test
    ).collectFile(storeDir: "${params.outdir}/assessment", keepHeader: true, skip: 1)
    
    // compute clustering metrics
    COMPUTE_CLUSTERING_METRICS(
        bins_fmt.bins,
        truth
    ).collectFile(storeDir: "${params.outdir}/assessment", keepHeader: true, skip: 1)

    if (params.with_checkv) {
        CHECKV_END_TO_END(
            bins_fmt.fasta.map{[[id: it[0].values().join('-')], it[1..-1]]},
            file(params.checkv.db, checkIfExists: true),
            [publish_dir: 'checkV']
        )
    }
}


workflow sim {
    fasta = Channel.fromPath("$params.root/*/assembly.fasta")
        .map{[[id: it.getParent().getName()], it]}
    coverage = Channel.fromPath("$params.root/*/coverage_contigs.h5")
        .map{[[id: it.getParent().getName()], it]}
    truth = file(params.truth, checkIfExists: false)
    
    binning(fasta, coverage, truth)
}

workflow single {
    // name = "${file(params.fasta).getParent().getName()}"
    name = 'Station-ALOHA'
    fasta = Channel.value([[id: "$name"], file(params.fasta, checkIfExists: true)])
    coverage =  Channel.value([[id: "$name"], file(params.coverage, checkIfExists: true)])
    truth = file(params.truth, checkIfExists: false)
    binning(fasta, coverage, truth)
}

workflow {
    params.coconet.args = "--n-train 10000"
    params.root = "${baseDir}/test"
    sim()
}
