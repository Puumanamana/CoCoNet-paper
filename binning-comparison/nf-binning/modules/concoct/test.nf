#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CONCOCT } from './process'


workflow test_concoct {
    def input = []
    input = [
        [ id: 'test' ],
        file("$baseDir/input/contigs.fa", checkIfExists: true),
        file("$baseDir/input/coverage.tsv", checkIfExists: true)
    ]

    CONCOCT (
        input,
        [ publish_dir:'test_concoct' ]
    )
}


workflow {
    test_concoct()
}
