#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {DOWNLOAD_VIRAL_REFSEQ;
         DOWNLOAD_TAXONOMY;
         GENERATE_CONFIG ;
         CAMISIM ;
         PROCESS_CAMISIM_OUTPUT ;
         SAMTOOLS_DEPTH;
         TO_H5} from './process'


// Main logic
workflow camisim {
    take:
    coverage_levels
    nb_samples
    nb_genomes
    nb_replicates

    main:
    db = DOWNLOAD_VIRAL_REFSEQ()
    tax_db = DOWNLOAD_TAXONOMY()

    configs = GENERATE_CONFIG(
        db,
        Channel.from(coverage_levels), 
        Channel.from(nb_samples),
        Channel.from(nb_genomes),
        Channel.from(1..nb_replicates)
    )

    simulation = CAMISIM(
        configs.camisim,
        tax_db
    )
    
    sim_info = PROCESS_CAMISIM_OUTPUT(
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
            .combine(configs.sizes, by: 0)
    )        
}
 

// Runner on paper settings
workflow all {
    camisim(
        [4, 10],
        [3, 15],
        [500, 2000],
        10
    )
}

workflow optimize {
    camisim(
        [4, 10],
        [3, 15],
        [500],
        1
    )
}

// Test runner
workflow test {
    camisim(
        [1],
        [2, 3],
        [5],
        2
    )
}
