process BAM_TO_ABD_TABLE{
    container 'quay.io/biocontainers/metabat2:2.15--h986a166_1'
    conda (params.conda ? 'bioconda::metabat2=2.15' : null)

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

    container 'nakor/coconet-paper-python'
    conda (params.conda ? 'pandas numpy h5py biopython' : null)

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

process MERGE_BINS {
    publishDir "${params.outdir}/merged_bins", mode: 'copy'
    container 'nakor/coconet-paper-python'
    conda (params.conda ? 'biopython pandas bioconda::blast' : null)
    
    input:
    tuple val(meta), path(bins), path(fasta)
    
    output:
    tuple val(meta), path("*-merged.fasta"), emit: fasta
    tuple val(meta), path("*-complete.csv"), emit: bins    
    script:
    """
    merge_bins.py \\
        --fasta $fasta \\
        --bins $bins \\
        --min-dtr-size 10 \\
        --max-dtr-size 300 \\
        --min-dtr-id 0.95
    """    
}

process COMPUTE_NN_METRICS {
    container 'nakor/coconet-paper-python'
    conda (params.conda ? 'scikit-learn pandas' : null)

    input:
    tuple val(meta), path(test)

    output:
    path "nn_scores.csv"

    script:
    """
    compute_nn_metrics.py --test $test --output nn_scores.csv --name $meta.id
    """
}

process COMPUTE_CLUSTERING_METRICS {
    tag {"${meta.id}"}
    container 'nakor/coconet-paper-python'
    conda (params.conda ? 'scikit-learn pandas' : null)

    input:
    tuple val(meta), path(bins)
    path truth

    output:
    path "clustering_scores.csv"

    script:
    truth_args = truth.isFile() ? "--truth ${truth}" : ""
    """
    compute_clustering_metrics.py $truth_args \\
        --bins ${bins.join(' ')} \\
        --output clustering_scores.csv \\
        --name $meta.id
    """
}
