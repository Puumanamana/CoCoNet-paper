process BAM_TO_ABD_TABLE{
    conda 'bioconda::metabat2=2.15'

    input:
    tuple val(meta), path(coverage)

    output:
    tuple val(meta), path("coverage_*.tsv")

    script:
    """
    jgi_summarize_bam_contig_depths \\
         --minMapQual ${params.min_map_qual} \\
         --percentIdentity ${params.min_map_id} \\
        -outputDepth coverage_${meta.id}.tsv 
    """
}

process REFORMAT_COVERAGE {
    publishDir "${params.outdir}/coverage_tables", mode: 'copy'
    conda 'pandas numpy h5py biopython'

    input:
    tuple val(meta), path(fasta), path(coverage)

    output:
    tuple val(meta), path("coverage_metabat2*.tsv"), emit: metabat2
    tuple val(meta), path("coverage_concoct*.tsv"), emit: concoct
    tuple val(meta), path("coverage_maxbin2*.tsv"), emit: maxbin2

    script:
    """
    h5_to_summary_table.py \\
        --abundance $coverage \\
        --fasta $fasta \\
        --suffix $meta.id
    """    
}


process BINS_TO_FASTA {
    publishDir "${params.outdir}/merged_bins", mode: 'copy'
    conda 'biopython'

    input:
    tuple val(meta), path(bins), path(fasta)
    
    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.csv"), emit: bins // corrected bins if singletons were set aside
    
    script:
    """
    bin_postprocessing.py \\
        --fasta $fasta \\
        --bins $bins
    """
}


process COMPUTE_SCORES {
    publishDir "${params.outdir}/scores", mode: 'copy'
    conda 'scikit-learn pandas'

    input:
    val(ds)
    path(bins)
    path truth

    output:
    tuple val(ds), path("*.csv")

    script:
    truth_args = truth.isFile() ? "" : "--truth ${truth}"
    """
    compute_metrics.py $truth_args \\
        --bins ${bins.join(' ')} \\
        --output scores-${ds}.csv
    """
}
