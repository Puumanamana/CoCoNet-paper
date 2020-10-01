include { initOptions; saveFiles ; getSoftwareName } from './functions'

process COCONET_PREPROCESS {
    tag {"${meta.id}"}
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container 'nakor/coconet'
    // Submit bioconda release once it's published
    // conda (params.conda ? "bioconda::coconet-binning=0.54" : null)
    
    input:
    tuple val(meta), file(fasta), file(coverage)
    val options

    output:
    tuple val(meta), path("output/*.fasta"), emit: fasta
    tuple val(meta), path("output/*.h5"), emit: h5
    tuple val(meta), path("output/*.log"), emit: log
    path "*.version.txt", emit: version

    script:
    def ioptions = initOptions(options)
    def software = getSoftwareName(task.process)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    def cov_list = coverage instanceof List ? coverage : [coverage]
    def cov_flag = cov_list[0].getExtension() == 'bam' ? "--bam" : "--h5"
    """
    coconet preprocess $ioptions.args \\
        --fasta $fasta \\
        $cov_flag ${cov_list.join(' ')} \\
        --output output \\
        --thread $task.cpus

    coconet --version > ${software}.version.txt
    """
}

