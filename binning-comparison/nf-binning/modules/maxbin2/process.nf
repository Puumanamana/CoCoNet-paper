include { initOptions; saveFiles ; getSoftwareName } from './functions'

/*
# Main parameters

-min_contig_length     minimum contig length. Default 1000
-prob_threshold        probability threshold for EM final classification. Default 0.9
*/

process MAXBIN2 {
    tag {"${meta.id}"}
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container 'nanozoo/maxbin2'
    // issue with the bioconda repository (missing lwp-simple perl library)
    conda (params.conda ? 'bioconda::maxbin2=2.2.7 perl-lwp-simple' : null)

    input:
    tuple val(meta), path(fasta), path(coverage)
    val options

    output:
    tuple val(meta), path('maxbin2-*.csv'), emit: bins
    path "*.version.txt", emit: version
    
    script:
    def ioptions = initOptions(options)
    def software = getSoftwareName(task.process)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def cov_arg = coverage.withIndex()
        .collect{v, i -> "-abund${i+1} ${v}"}
        .join(' ')
        .replaceFirst('abund1', 'abund')
    """
    run_MaxBin.pl $ioptions.args \\
        -contig $fasta $cov_arg \\
        -out maxbin-$prefix \\
        -thread $task.cpus

    # Make assignment file
    for fa in `ls maxbin-${prefix}*.fasta`; do
        bin_id=\$(echo \$fa | cut -d '.' -f2 | sed 's/^0*//')
        grep '^>' \$fa | cut -c 2- | awk -v b=\$bin_id '{print \$1","b}' >> maxbin2-${prefix}.csv
    done

    run_MaxBin.pl -v | head -1 > ${software}.version.txt
    """
}
