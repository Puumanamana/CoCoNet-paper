process GENERATE_CONFIG {
    tag "$meta.id"
    label 'low_computation'
    conda 'pandas biopython configparser'
	
    input:
    each xcoverage
    each n_samples
    each n_genomes
	each replicate

    output:
    tuple val(sim_params), file('*.{ini,tsv}')
    
    script:
    meta = [
        id: "camisim_lmu.${params.log_mu}-lsig.${params.log_sigma}-ns.${n_samples}-cov.${xcoverage}X-ng.${n_genomes}-r.${replicate}",
        coverage: xcoverage,
        n_samples: n_samples,
        n_genomes: n_genomes,
        replicate: replicate,
    ]
    """
    generate_camisim_metadata.py \
        --name $meta.id
        --cov_lvl $xcoverage \
        --n_samples $n_samples \
        --n_contigs $n_genomes \
        --cami-data $params.cami_data \
        --log-mu 1 \
        --log-sigma 2
    """
}

process CAMISIM {
    tag "$meta.id"
    label 'high_computation'
    container 'nakor/coconet-paper-camisim'
    
    input:
    tuple val(meta), path(config_files)
    path contig_db
    
    output:
    set val(meta), file("camisim_*/gsa_pooled.fasta"), emit: assembly
    set val(meta), file("camisim_*/*sample*/bam/*.bam"), emit: bam

	script:
    """
    metagenomesimulation.py -s0 -p $task.cpus config.ini
    """
}

process GENERATE_METADATA {
    publishDir "${params.outdir}/${id}", mode: "copy"
    label 'low_computation'
    conda 'pandas'

	input:
    tuple val(id), path(fasta)

	output:
    tuple val(id), path('metadata.csv')
    path 'assembly.fasta'

	script:
    """
    import pandas as pd

    ctg_info = []
    for line in open("$fasta"):
        if line.startswith('>'):
            (name, start, end, size) = [x[i] for x in line.strip().split('_') for i in range(0, 7, 2)]
            ctg_info.append([name, start, end, size]

    ctg_info = pd.DataFrame(ctg_info, columns=['V_id', 'start', 'end', 'size'])

    suffixes = ctg_info.groupby('V_id').cumcount().astype(str)
    ctg_info['C_id'] = ctg_info.V_id + '|' + suffixes

    ctg_info.astype(int).to_csv('metadata.csv', index=False)

    i = 0
    writer = open('assembly.fasta', 'w')
    for line in open(assembly_file):
         if not line.startswith('>'):
             writer.write(line)
         else:
             writer.write(ctg_info.C_id[i])
             i += 1
    writer.close()
    """
}

process SAMTOOLS_DEPTH {
    label 'medium_computation'
    tag "${meta.id}_${genome}"
    conda 'samtools'

	input:
    tuple val(meta), val(genome), path(bams)

	output:
    tuple val(meta), val(genome), path("*.txt")

    script:
    """
    samtools depth ${bams.join(' ')} > ${contig}.txt
    """
}

process TO_H5 {
    tag {"$meta.id"}
    publishDir "${params.outdir}/${id}", mode: "copy"
    label 'medium_computation'
    conda 'pandas h5py'
	
    input:
    tuple val(meta), file(depth), file(meta_file)

    output:
    file("coverage_virus.h5")
    file("coverage_contigs.h5")

	script:
    """
    #!/usr/bin/env python

    from pathlib import Path
    import pandas as pd
    import h5py

    metadata = pd.read_csv("${meta_file}")

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



