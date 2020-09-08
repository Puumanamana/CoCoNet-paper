nextflow.enable.dsl=2

include { REFORMAT_COVERAGE; BINS_TO_FASTA } from './util/process'
include { COCONET_PREPROCESS } from './modules/binning/coconet/preprocess/process'
include { COCONET_RUN } from './modules/binning/coconet/run/process' 
include { CONCOCT } from './modules/binning/concoct/process'
include { METABAT2 } from './modules/binning/metabat2/process'
include { MAXBIN2 } from './modules/binning/maxbin2/process'


workflow binning {
    take:
    folder

    main:
    ds = "${folder.getName()}"
    fasta = Channel.fromPath("${folder}/*.fasta").map{[ [id: ds], it]}
    bams = Channel.fromPath("${folder}/*.bam").map{[ [id: ds], it]}
    h5 = Channel.fromPath("${folder}/*.{h5,hdf5}").map{[ [id: ds], it]}

    preproc_args = ["--min-ctg-len", "$params.min_ctg_len",
                    "--min-mapping-quality", "$params.min_mapq",
                    "--min-aln-coverage", "$params.min_map_id",
                    "--fl-range", "$params.fl_range"]

    coconet_args = preproc_args + ["--min-prevalence", "${params.coconet.min_prevalence}",
                                   "--theta", "${params.coconet.theta}",
                                   "--gamma1", "${params.coconet.gamma1}",
                                   "--gamma2", "${params.coconet.gamma2}",
                                   "--features", "${params.coconet.features}"]

    // inputs preprocessing
    inputs = COCONET_PREPROCESS(
        fasta.join(bams.mix(h5)),
        [publish_dir: 'preprocessing',
         args: preproc_args.join(' ')]
    )
    composition = inputs.fasta
    coverage = REFORMAT_COVERAGE(inputs.h5)

    // binning
    coconet_bins = COCONET_RUN(
        composition.join(inputs.h5),
        [publish_dir: 'CoCoNet',
         args: coconet_args.join(' ')]
    ).bins
    metabat2_bins = METABAT2(
        composition.join(coverage.metabat2),
        [publish_dir: 'Metabat2',
         args: "--minContig ${params.min_ctg_len}"]
    ).bins
    concoct_bins = CONCOCT(
        composition.join(coverage.concoct),
        [publish_dir: 'CONCOCT',
         args: "--clusters ${params.concoct.max_clusters} --length_threshold ${params.min_ctg_len}"]
    ).bins
    maxbin2_bins = MAXBIN2(
        composition.join(coverage.maxbin2),
        [publish_dir: 'MaxBin2',
         args: "-min_contig_length ${params.min_ctg_len}"]
    ).bins

    // reformating
    BINS_TO_FASTA(
        coconet_bins.mix(metabat2_bins, concoct_bins, maxbin2_bins).join(fasta)
    )
}


workflow all_sims {
    Channel.fromPath("${params.root_dir}/*", type: 'dir') | binning
}

workflow {
    binning(file("${baseDir}/test"))
}
