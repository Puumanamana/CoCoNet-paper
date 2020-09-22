#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {GENERATE_CONFIG ; CAMISIM ; GENERATE_METADATA; SAMTOOLS_DEPTH; TO_H5} from './process'

workflow {
    db = file(params.db, checkIfExists: true)
    
    configs = GENERATE_CONFIG(
        db,
        Channel.from(4, 10), // coverage
        Channel.from(5, 20), // number of samples
        Channel.from(500, 2000), // number of genomes
        Channel.from(0..10) // number of replicates
    )

    simulation = CAMISIM(
        configs,
        db
    )
    
    sim_info = GENERATE_METADATA(
        simulation.assembly
    )

    coverage_txt = SAMTOOLS_DEPTH(
        simulation.bam
            .transpose()
            .map{[it[0], it[1].toString().find(/V\d+/), it[1]]}
            .groupTuple(by: [0, 1])
    )

    TO_H5(
        coverage_txt
            .groupTuple(by: 0)
            .collect{it[1][1]}
            .combine(sim_info, by: 0)
    )
    
}

