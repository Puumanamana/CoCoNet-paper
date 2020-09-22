#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {DOWNLOAD_VIRAL_REFSEQ; DOWNLOAD_TAXONOMY;
         GENERATE_CONFIG ; CAMISIM ;
         GENERATE_METADATA;
         SAMTOOLS_DEPTH; TO_H5} from './process'

workflow camisim {
    take:
    xcov_lvls
    nb_samples
    nb_genomes
    nb_replicates

    main:
    db = DOWNLOAD_VIRAL_REFSEQ()
    tax_db = DOWNLOAD_TAXONOMY()
    
    configs = GENERATE_CONFIG(
        db,
        Channel.from(xcov_lvls), // coverage
        Channel.from(nb_samples), // number of samples
        Channel.from(nb_genomes), // number of genomes
        Channel.from(0..nb_replicates) // number of replicates
    )

    simulation = CAMISIM(
        configs,
        tax_db
    )
    
    sim_info = GENERATE_METADATA(
        simulation.assembly
    )

    coverage_txt = SAMTOOLS_DEPTH(
        simulation.bam
            .transpose()
            .map{[it[0], it[1].getSimpleName().tokenize('-')[0], it[1]]}
            .groupTuple(by: [0, 1])
    )
    
    TO_H5(
        coverage_txt
            .groupTuple(by: 0)
            .combine(sim_info.table, by: 0)
    )        
}
 
workflow all {
    camisim(
        [4, 10],
        [5, 20],
        [500, 2000],
        10
    )
}

workflow test {
    camisim(
        [1],
        [2, 3],
        [5],
        2
    )
}
