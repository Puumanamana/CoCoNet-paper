#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAXBIN2 } from './process'


workflow test_maxbin2 {
    def input = []
    input = [
        [ id: 'test' ],
        file("$baseDir/input/contigs.fa", checkIfExists: true),
        file("$baseDir/input/coverage.tsv", checkIfExists: true)
    ]

    MAXBIN2 (
        input,
        [ publish_dir:'test_maxbin2' ]
    )
}

workflow {
    test_maxbin2()
}
