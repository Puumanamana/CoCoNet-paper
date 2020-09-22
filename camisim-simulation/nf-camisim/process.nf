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
    label 'low_computation'
    
    output:
    path '*.fna'

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{1,2}.1.genomic.fna.gz
    zcat *.fna.gz | sed '/^>/! s/N//g' > viral_refseq.fna
    """
}

process GENERATE_CONFIG {
    tag "$meta.id"

    container 'nakor/coconet-paper-python'
    label 'low_computation'

    input:
    path db
    each xcoverage
    each n_samples
    each n_genomes
    each replicate

    output:
    tuple val(meta), path('*.{ini,tsv}'), path('source-genomes')

    script:
    meta = [
        id: "camisim.lmu-${params.log_mu}.lsig-${params.log_sigma}.ns-${n_samples}.cov-${xcoverage}X.ng-${n_genomes}.${replicate}",
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
        --cami-data $params.cami_data \\
        --log-mu 1 \\
        --log-sigma 2
    """
}

process CAMISIM {
    tag "$meta.id"
    publishDir "$params.outdir/$meta.id/bam"

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
    for bam in `ls camisim*/*sample*/bam/*.bam` ; do
        [[ "\$bam" =~ sample_[0-9]* ]] && sample_nb=\$BASH_REMATCH || echo "no sample id found" ||  exit 1
        mv \$bam \${bam%.*}-\${sample_nb}.bam
    done
    """
}

process GENERATE_METADATA {
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
    #!/usr/bin/env python

    import pandas as pd
    import re

    pattern = re.compile(r'>(.*)_from_([0-9]+)_to_([0-9]+)_total_([0-9]+)')
    ctg_info = []
    for line in open("$fasta"):
        if line.startswith('>'):
            (name, start, end, size) = re.findall(pattern, line)[0]
            ctg_info.append([name, int(start), int(end), int(size)])

    ctg_info = pd.DataFrame(ctg_info, columns=['V_id', 'start', 'end', 'size'])

    suffixes = ctg_info.groupby('V_id').cumcount().astype(str)
    ctg_info['C_id'] = ctg_info.V_id + '|' + suffixes

    ctg_info.to_csv('metadata.csv', index=False)

    i = 0
    writer = open('assembly.fasta', 'w')
    for line in open("$fasta"):
         if not line.startswith('>'):
             writer.write(line)
         else:
             writer.write(ctg_info.C_id[i])
             i += 1
    writer.close()
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
    tuple val(meta), val(genome), path("*.txt")

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
    tuple val(meta), path(depth), path(sim_info)

    output:
    path("coverage_virus.h5")
    path("coverage_contigs.h5")

    script:
    """
    #!/usr/bin/env python

    from pathlib import Path
    import pandas as pd
    import h5py

    metadata = pd.read_csv("$sim_info")

    info = metadata.groupby("V_id").agg(list)

    cov_vir_h5 = h5py.File("coverage_virus.h5","w")
    cov_ctg_h5 = h5py.File("coverage_contigs.h5","w")

    for filename in Path('.').glob('*.txt'):
        virus = filename.stem

        if virus not in info.index:
            continue

        coverage = (
            pd.read_csv(filename, usecols=lambda x: x!=0, dtype=int, index_col=0)
            .set_axis(coverage.index - 1, axis=0) # since samtools positions are 1-based
            .set_axis([f'sample_{i+1}' for i in range($meta.n_samples)], axis=1)
            .reindex(index=range(metadata.loc[virus, 'size']))
            .fillna(0)
            .to_numpy()
        )

        cov_vir_h5.create_dataset(virus, data=coverage)

        for ctg, start, end in zip(*info.loc[virus,["C_id","start","end"]]):
            cov_ctg_h5.create_dataset(ctg, data=matrix[:,start:end+1])

    cov_ctg_h5.close()
    cov_vir_h5.close()
    """
}
