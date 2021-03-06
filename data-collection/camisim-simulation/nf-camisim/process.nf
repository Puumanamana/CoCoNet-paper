process DOWNLOAD_TAXONOMY {
    label 'low_computation'
    
    output:
    path '*.tar.gz'

    script:
    """
    mkdir NCBI
    wget -qO- ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz \
    | tar xz -C NCBI 
    tar czf new_taxdump.tar.gz NCBI
    """
}

process DOWNLOAD_VIRAL_REFSEQ {
    publishDir "$params.outdir/refseq-db"
    label 'low_computation'
    container 'nakor/coconet-paper-python'
    
    output:
    path '*ACGT.fna'

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{1,2}.1.genomic.fna.gz
    zcat *.fna.gz | sed '/^>/! s/N//g' > viral-refseq.fna
    clean_db.py --fasta viral-refseq.fna --output viral-refseq-ACGT.fna
    """
}

process GENERATE_CONFIG {
    tag "$meta.id"

    container 'nakor/coconet-paper-python'
    label 'low_computation'

    input:
    path db
    each n_genomes
    each n_samples
    each xcoverage
    each replicate

    output:
    tuple val(meta), path('*.{ini,tsv}'), path('source-genomes'), emit: camisim
    tuple val(meta), path('genome_sizes.csv'), emit: sizes

    script:
    meta = [
        id: "camisim_${n_genomes}-genomes_${n_samples}-samples_${xcoverage}X_${replicate}",
        coverage: xcoverage,
        n_samples: n_samples,
        n_genomes: n_genomes,
        replicate: replicate,
    ]
    """
    generate_camisim_metadata.py \\
        --db $db \\
        --min-genome-size $params.min_genome_size \\
        --name $meta.id \\
        --cov_lvl $xcoverage \\
        --n_samples $n_samples \\
        --n_genomes $n_genomes \\
        --cami-data ${params.camisim.data} \\
        --log-mu ${params.camisim.log_mu} \\
        --log-sigma ${params.camisim.log_sigma}
    """
}

process CAMISIM {
    tag "$meta.id"
    publishDir "$params.outdir/$meta.id"

    container 'nakor/coconet-paper-camisim'
    label 'high_computation'

    input:
    tuple val(meta), path(config_files), path(source_genomes)
    path taxdump

    output:
    tuple val(meta), path("camisim*/gsa_pooled.fasta"), emit: assembly
    tuple val(meta), path("camisim*/*sample*/bam/*.bam"), emit: bam

    script:
    """
    #!/usr/bin/env bash

    metagenomesimulation.py -s0 -p $task.cpus config.ini

    # Rename bam files
    rename-camisim.py
    """
}

process PROCESS_CAMISIM_OUTPUT {
    tag "$meta.id"
    publishDir "$params.outdir/$meta.id", mode: "copy"

    container 'nakor/coconet-paper-python'
    label 'low_computation'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('metadata.csv'), emit: table
    path 'assembly.fasta', emit: assembly

    script:
    """
    process_camisim_output.py --fasta $fasta
    """
}

process SAMTOOLS_DEPTH {
    tag "${meta.id}_${genome}"
    publishDir "$params.outdir/$meta.id/txt"

    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'
    label 'low_computation'

    input:
    tuple val(meta), val(genome), path(bams)

    output:
    tuple val(meta), path("*.txt")

    script:
    """
    samtools depth ${bams.join(' ')} > ${genome}.txt
    """
}

process TO_H5 {
    tag {"$meta.id"}
    publishDir "$params.outdir/$meta.id", mode: "copy"

    label 'medium_computation'
    container 'nakor/coconet-paper-python'

    input:
    tuple val(meta), path(depth), path(sim_info), path(genome_sizes)

    output:
    path("coverage_virus.h5")
    path("coverage_contigs.h5")

    script:
    """
    depth_to_h5.py --metadata $sim_info --genome-sizes $genome_sizes --n-samples $meta.n_samples
    """
}
