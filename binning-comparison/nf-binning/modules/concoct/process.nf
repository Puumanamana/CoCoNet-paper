include { initOptions; saveFiles ; getSoftwareName } from './functions'

/*
Main parameters:
 --length_threshold: 
 --clusters: specify maximal number of clusters for VGMM, default 400
 --no_original_data: By default the original data is saved to disk. 
                     For big datasets, especially when a large k is used 
                     for compositional data, this file can become very large. 
                     Use this tag if you don't want to save the original data.
*/

process CONCOCT {
    tag {"${meta.id}"}
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
    
    container "quay.io/biocontainers/concoct:1.1.0--py38h7be5676_2"
    conda (params.conda ? "bioconda::concoct=1.1.0" : null)

    input:
    tuple val(meta), path(fasta), path(coverage)
    val options

    output:
    tuple val(meta), path("concoct_${meta.id}.csv"), emit: bins
    path "*.version.txt", emit: version

    script:
    def ioptions = initOptions(options)
    def software = getSoftwareName(task.process)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"    
    """
    concoct $ioptions.args \\
        --threads $task.cpus \\
        --composition_file ${fasta} \\
        --coverage_file ${coverage}

    tail -n+2 clustering_gt*.csv > concoct_${meta.id}.csv

    concoct --version > ${software}.version.txt
    """
}
