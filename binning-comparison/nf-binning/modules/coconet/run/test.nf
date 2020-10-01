#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COCONET_RUN } from './process'

workflow test_coconet_run {
    def input = []
    input = [
        [ id: 'test' ],
        file("$baseDir/input/contigs.fa", checkIfExists: true),
        file("$baseDir/input/sample*.bam", checkIfExists: true)
    ]

    COCONET_RUN (
        input,
        [ publish_dir:'test_coconet_run' ]
    )
}


workflow {
    test_coconet_run()
}
