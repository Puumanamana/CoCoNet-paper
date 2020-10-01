#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METABAT2 } from './process'


workflow test_metabat2 {
    def input = []
    input = [
        [ id: 'test' ],
        file("$baseDir/input/contigs.fa", checkIfExists: true),
        file("$baseDir/input/depth.txt", checkIfExists: true)
    ]

    METABAT2 ( input,
              [ publish_dir:'test_metabat2' ] )
}


workflow {
    test_metabat2()
}
