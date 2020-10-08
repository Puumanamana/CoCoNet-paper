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
    n_genomes
    n_samples
    coverage_levels
    n_replicates

    main:
    db = DOWNLOAD_VIRAL_REFSEQ()
    tax_db = DOWNLOAD_TAXONOMY()

    configs = GENERATE_CONFIG(
        db, n_genomes, n_samples, coverage_levels, 1..n_replicates
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
 

// Runner
// For testing, use default parameters
workflow {
    camisim(
        params.n_genomes.tokenize(','),
        params.n_samples.tokenize(','),
        params.x_coverage.tokenize(','),
        params.n_replicates
    )
}
