#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


GENOME_LENGTHS = Channel
    .fromPath(params.db)
    .splitFasta( record: [id: true, seqString: true ])
    .map{ it -> tuple( it.id, it.seqString.length() ) }

process GENERATE_CONFIG {
    label 'low_computation'
	
    input:
    each coverage_level
    each n_samples
    each n_genomes
	each replicate

    output:
    tuple val(coverage_level), val(n_samples), val(n_genomes), val(replicate), file('config.ini'), file("id_to_genome.tsv"), file('metadata_camisim.tsv')
    
    script:
    """
    generate_camisim_metadata.py \
             --cov_lvl ${coverage_level} \
             --n_samples ${n_samples} \
             --n_contigs ${n_genomes} \
             --db ${params.db} \
             --cami_src ${params.camisim_path}
    """
}

process CamisimRunner {
    tag "Camisim_ng=${n_genomes}_cov=${coverage_level}_nsamples=${n_samples}"
    label 'high_computation'	
    label 'python2'
    
    input:
    set coverage_level, n_samples, n_genomes, replicate, file(cfg), file(idtg), file(meta)
    file(contigs)
    
    output:
    set val(meta), file("sim_ng*/gsa_pooled.fasta")
    set val(meta), file("sim_ng*/*sample*/bam/*.bam")

	script:
    """
    metagenomesimulation.py -s0 -p ${task.cpus} config.ini
    """
}

process GenerateMetadataFile {
    publishDir "${params.outdir}/${id}", mode: "copy"
    label 'low_computation'
    label 'python3'

	input:
    set id, file(fasta) from CAMISIM_ASSEMBLY

	output:
    set id, file('metadata.csv') into METADATA_FILES
    file('assembly.fasta')

	script:
    """
    # id = ${id}
    python ${workflow.projectDir}/extract_contig_info.py ${fasta}
    """
}

// CAMISIM_COVERAGE emits items of the shape [id, xxx/sample_i/V_j.bam]
// Reshapes the channel with items of the form [ id, sample, virus, file ]
CAMISIM_COVERAGE_SPLIT = CAMISIM_COVERAGE
    .map{[it[0],
		  it[1].collect { [it.toString().find(/sample_\d+/),
						   it.toString().find(/V\d+/),
						   it ]}]}
    .transpose()
    .map{ [it[0], it[1][0], it[1][1], it[1][2]]}

process Samtools {
    label 'low_computation'
    label 'samtools'
    tag "samtools_${id}"

	input:
    set id, sample, contig, file(bam) from CAMISIM_COVERAGE_SPLIT

	output:
    set contig, id, file("*.txt") into DEPTH_TXT

    script:
    """
    #!/usr/bin/env bash

    # id = ${id}
    samtools depth $bam > "${sample}_${contig}.txt"
    """
}

DEPTH_TXT_RESHAPED = DEPTH_TXT
	.groupTuple(by: [0,1])

process TxtToNpy {
    tag "${id}_${contig}"
    label 'medium_computation'
    label 'python3'

	input:
    set val(contig), val(id), file(f), val(genome_len) from DEPTH_TXT_RESHAPED.combine(GENOME_LENGTHS,by: 0)

	output:
    set val(id), file("*.npy") into DEPTH_NPY

    script:
    """
    #!/usr/bin/env python
    # id=${id}

    import pandas as pd
    import numpy as np

    files = sorted([ "${f.join('","')}" ])

    datasets = [ pd.read_csv(f, header=None, sep='\\t', names=["acc","pos","depth"]) 
                 for f in sorted(files) ]

    coverage = np.zeros([len(files), ${genome_len}],dtype=np.uint32)

    for i,dataset in enumerate(datasets):
        if dataset.size > 0:
            coverage[i,dataset.pos.values-1] = dataset.depth.values

    np.save("${contig}.npy",coverage)
    """
}

DEPTH_NPY_GROUPED = DEPTH_NPY.groupTuple()

process NpyToH5 {
    tag "${id}"
    publishDir "${params.outdir}/${id}", mode: "copy"
    label 'high_computation'
    label 'python3'
	
    input:
    set id, file(f), file(meta_file) from DEPTH_NPY_GROUPED.join(METADATA_FILES)

    output:
    file("coverage_virus.h5")
    file("coverage_contigs.h5")

	script:
    """
    #!/usr/bin/env python

    id = ${id}
    from glob import glob
    import numpy as np
    import pandas as pd
    import h5py

    metadata = pd.read_csv("${meta_file}")

    info = metadata.groupby("V_id").agg(list)

    cov_vir_h5 = h5py.File("coverage_virus.h5","w")
    cov_ctg_h5 = h5py.File("coverage_contigs.h5","w")
    
    for filename in glob("*.npy"):
        virus = filename.split(".")[0]
        matrix = np.load(filename)

        if virus not in info.index:
            continue
        cov_vir_h5.create_dataset(virus,data=matrix)

        for ctg,start,end in zip(*info.loc[virus,["C_id","start","end"]]):
            cov_ctg_h5.create_dataset(ctg,data=matrix[:,start:end+1])

    cov_ctg_h5.close()
    cov_vir_h5.close()
    """
}



