nextflow.enable.dsl = 2

include { DOWNLOAD_ILLUMINA_READS; DOWNLOAD_ONT_ASSEMBLY;
         TRIMMING_FASTP; COASSEMBLY_METASPADES;
         BWA_INDEX; BWA_MEM;
         DEDUPLICATION_ONT; MAP_ILLUMINA_ONT} from "./process"

workflow {

    // Generate assembly from illumina reads
    accessions = Channel.from(["SRR5677468", "SRR8811962", "SRR8811963"])

    reads = DOWNLOAD_ILLUMINA_READS(accessions)
    trimmed_reads = TRIMMING_FASTP(reads)
    illumina_assembly = COASSEMBLY_METASPADES(trimmed_reads.map{it[1]}.collect())

    // Compute coverage on assembly
    ref = BWA_INDEX(illumina_assembly)
    coverage = BWA_MEM(ref, trimmed_reads)
    
    // Download and merge ONT assembly from NCBI
    // draft_genomes = DOWNLOAD_ONT_ASSEMBLY("PRJNA529454")
    // draft_genomes_unique = DEDUPLICATION_ONT(draft_genomes)

    // Map illumina assembly on ONT assembly
    // MAP_ILLUMINA_ONT(illumina_assembly, draft_genomes_unique)
    
}
