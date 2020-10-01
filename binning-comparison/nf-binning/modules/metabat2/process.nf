include { initOptions; saveFiles ; getSoftwareName } from './functions'

/*
# Main parameters

-m [ --minContig ] arg (=2500)    Minimum size of a contig for binning (should be >=1500).
-s [ --minClsSize ] arg (=200000) Minimum size of a bin as the output.
-l [ --onlyLabel ]                Output only sequence labels as a list in a column without sequences.
--saveCls                         Save cluster memberships as a matrix format
--unbinned                        Generate [outFile].unbinned.fa file for unbinned contigs
--noBinOut                        No bin output. Usually combined with --saveCls to check only contig memberships
--seed arg (=0)                   For exact reproducibility. (0: use random seed)
*/

process METABAT2 {
    tag {"${meta.id}"}
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
    
    container 'quay.io/biocontainers/metabat2:2.15--h986a166_1'
    conda (params.conda ? "bioconda::metabat2=2.15" : null)

    input:
    tuple val(meta), path(fasta), path(coverage)
    val options

    output:
    tuple val(meta), path("metabat2-${meta.id}.csv"), emit: bins
    path "*.version.txt", emit: version
    
    script:
    def ioptions = initOptions(options)
    def software = getSoftwareName(task.process)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    """
    gzip -c $fasta > contigs.fasta.gz && \\
    metabat2 $ioptions.args \\
        --inFile contigs.fasta.gz \\
        --abdFile $coverage \\
        --outFile metabat2-${meta.id} \\
        --unbinned \\
        --numThreads $task.cpus

    # Make assignment file
    for f in `ls *.fa`; do
        bin=`echo \$f | cut -d'.' -f2`

        if [[ "\$bin" =~ ^[0-9]+\$ ]]; then
            grep '^>' \$f | cut -c2- \\
                | awk -v b=\$bin '{print \$1","b}' \\
                >> "metabat2-${meta.id}.csv"
        else
            grep '^>' \$f | cut -c2- \\
                | awk -v b=\$bin '{print \$1","b"-"NR}' \\
                >> "metabat2-${meta.id}.csv"
        fi
    done

    rm contigs.fasta.gz

    metabat2 -h 2>&1 | head | grep version > ${software}.version.txt
    """
}

