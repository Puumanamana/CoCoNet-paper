process DOWNLOAD_ILLUMINA_READS {
    tag {"$sample"}
    publishDir "$params.outdir/illumina-reads/raw", mode: "copy"
    label "medium_computation"

    input:
    val sample
    
    output:
    tuple val(sample), path("*_{1,2}.fastq.gz")

    script:
    """
    # fasterq-dump fails for some accesions (github issue #318)
    # To circumvent this issue, we can use prefetch first

    prefetch $sample \\
    && fasterq-dump --threads $task.cpus --mem ${task.memory.getGiga()}GB \\
        --outdir . \\
        --split-files \\
        ./$sample

    pigz -p $task.cpus ${sample}_1.fastq
    pigz -p $task.cpus ${sample}_2.fastq
    """
}


process DOWNLOAD_ONT_ASSEMBLY {
    publishDir "$params.outdir/ONT-assemblies", mode: "copy"

    input:
    val accession
    
    output:
    path "*.fasta"

    script:
    """
    esearch -db nuccore -query "$accession" \\
        | efetch -format fasta \\
        > assemblies_ONT-concat.fasta
    """
}


process TRIMMING_FASTP {
    tag {"$sample"}
    publishDir "$params.outdir/illumina-reads/trimmed", mode: "copy"

    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("*_trimmed_{1,2,unpaired}.fastq.gz")

    script:
    """
    fastp \\
    --trim_poly_g \\
    --average_qual 25 \\
    --length_required 20 \\
    --cut_front --cut_tail \\
    --cut_mean_quality 25 --cut_window_size 10 \\
    --trim_front1 5 --trim_front2 5 \\
    --trim_tail1 5 --trim_tail2 5 \\
    -i ${fastqs[0]} \\
    -I ${fastqs[1]} \\
    -o ${sample}_trimmed_1.fastq.gz \\
    -O ${sample}_trimmed_2.fastq.gz \\
    --unpaired1 ${sample}_trimmed_unpaired.fastq.gz    
    """
}

process COASSEMBLY_METASPADES {
    publishDir "$params.outdir/assembly", mode: "copy"
    label "high_computation"
    
    input:
    path reads

    output:
    path "output/scaffolds.fasta"

    script:
    """
    cat *_1.fastq.gz > all_1.fastq.gz
    cat *_2.fastq.gz > all_2.fastq.gz
    cat *_unpaired.fastq.gz > all_unpaired.fastq.gz

    spades.py --meta -1 all_1.fastq.gz -2 all_2.fastq.gz -s all_unpaired.fastq.gz -o output

    rm -f all_*.fastq.gz
    """
}


process BWA_INDEX {

    input:
    path assembly

    output:
    path "*"

    script:
    """
    bwa index -p reference $assembly
    """
}

process BWA_MEM {
    tag {"$sample"}
    publishDir "$params.outdir/coverage", mode: "copy"
    label "medium_computation"
    
    input:
    path index
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.bam"), path("*.bai")

    script:
    """
    bwa mem -t $task.cpus -a -M reference $reads \\
        | samtools view -bh -@ $task.cpus \\
        | samtools sort -@ $task.cpus -o ${sample}.bam

    samtools index -@ $task.cpus ${sample}.bam ${sample}.bai
    """
}

process DEDUPLICATION_ONT {
    publishDir "$params.outdir/ONT-assemblies", mode: "copy"
    label "medium_computation"

    input:
    path assembly

    output:
    path "*.fna"

    script:
    """
    Cluster_genomes.pl -f $assembly -c 80 -i 95
    """
}

process MAP_ILLUMINA_ONT {
    publishDir "$params.outdir/mapping_illumina-vs-ONT", mode: "copy"
    label "medium_computation"

    input:
    path illumina_assembly
    path ont_assembly

    output:
    path "alignments.tsv", emit: minimap2
    path "mapping.csv", emit: bins

    script:
    """
    # 1) minimap2 alignment --> output in PAF format
    # 2) Filter alignments:
    # - mapq > 30 (column #12)
    # - query coverage > 95% (column #10/#2)
    # - no secondary alignments (flag tp:A:S)

    minimap2 $ont_assembly $illumina_assembly -t $task.cpus \\
        | grep -v 'tp:A:S'
        | awk '\$12 > 30 && \$10/\$2 > 0.95'
        > alignments.tsv

    # 3) Make mapping file
    cut -f1,6 alignments.tsv | sed 's/\\t/,/g' > mapping.csv

    # 4) Check if some illumina contigs map on multiple ONT contigs (case not handled)
    cut -f1 alignments.tsv | sort | uniq -c | awk '\$1 > 1' > multimappers.txt
    [ -s multimappers.txt ] && echo "Some contigs are mapping on multiple references. Aborting." && exit 1
    """
}
