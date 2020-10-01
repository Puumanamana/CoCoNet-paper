#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COCONET_PREPROCESS } from './process'

workflow test_coconet_preprocess {
    def input = []
    input = [
        [ id: 'test' ],
        file("$baseDir/input/contigs.fa", checkIfExists: true),
        file("$baseDir/input/sample*.bam", checkIfExists: true)
    ]

    COCONET_PREPROCESS (
        input,
        [ publish_dir:'test_coconet_preprocess' ]
    )
}


workflow {
    test_coconet_preprocess()
}
