nextflow.enable.dsl = 2

include { DOWNLOAD_ILLUMINA_READS; DOWNLOAD_ONT_ASSEMBLY;
         TRIMMING_FASTP; COASSEMBLY_METASPADES;
         BWA_INDEX; BWA_MEM;
         DEDUPLICATION_ONT; MAP_ILLUMINA_ONT} from "./process"

workflow illumina {
    take:
    accessions

    main:
    // Download
    reads = DOWNLOAD_ILLUMINA_READS(accessions)
    // Filter and assemble
    trimmed_reads = TRIMMING_FASTP(reads)
    illumina_assembly = COASSEMBLY_METASPADES(trimmed_reads.map{it[1]}.collect())
    // Compute coverage on assembly
    ref = BWA_INDEX(illumina_assembly)
    coverage = BWA_MEM(ref, trimmed_reads)

    emit:
    assembly = illumina_assembly
    coverage = coverage.bam_bai
}

workflow ont {
    take:
    project_accession
    illumina_assembly

    main:
    draft_genomes = DOWNLOAD_ONT_ASSEMBLY(project_accession)
    draft_genomes_unique = DEDUPLICATION_ONT(draft_genomes)
    MAP_ILLUMINA_ONT(illumina_assembly, ont_assembly)
}

workflow {
    // Generate assembly from illumina reads

    if (params.run_assembly) {
        accessions = Channel.from(["SRR5677468", "SRR8811962", "SRR8811963"])
        illumina_assembly = illumina(accessions).assembly
    } else {
        illumina_assembly = file(params.assembly, checkIfExists: true)
    }

    ont("PRJNA529454", illumina_assembly)
}
