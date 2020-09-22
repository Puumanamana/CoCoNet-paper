#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {GENERATE_CONFIG ; CAMISIM ; GENERATE_METADATA; SAMTOOLS_DEPTH; TO_H5} from process

workflow {
    configs = GENERATE_CONFIG(
        (4, 10), // coverage
        (5, 20), // number of samples
        (500, 2000), // number of genomes
        (0..10) // number of replicates
    )

    simulation = CAMISIM(
        configs,
        file(params.db, checkIfExists: true)
    )
    
    GENERATE_METADATA(
        simulation.assembly
    )

    coverage_txt = SAMTOOLS_DEPTH(
        simulation.bam
            .transpose()
            .map{[it[0], it[1].toString().find(/V\d+/), it[1]]}
            .groupTuple(by: [0, 1])
    )

    TO_H5(
        coverage_txt.groupTuple(by: 0).collect{it[1][1]}
    )
    
}

