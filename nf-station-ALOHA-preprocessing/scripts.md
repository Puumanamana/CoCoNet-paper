# Before starting

Global variables:

```bash
export SAMPLES="SRR5677468 SRR8811962 SRR8811963"
export THREADS=10
export MEM=500G
```


# Step 1: Download the illumina reads

Accessions:
- SRR5677468 (15m)
- SRR8811962 (117m)
- SRR8811963 (250m)

Conda environment: `conda install -c bioconda sra-tools=2.10.8`

Commands:
```bash
mkdir raw
fasterq-dump --threads $THREAD --mem $MEM --split-files $SAMPLES --outdir raw
for f in $(ls raw/*.fastq); do
    gzip $f
done
```

# Step 2: Trimming (fastp)

Parameters:
- Crop reads at beginning and end: 5 bp
- Trimming with sliding window of size 10 when quality < 25
- Minimum read length after trimming: 20 bp
- Trim poly-G sequences
    
Conda environment: `conda install -c bioconda -c conda-forge parallel fastp=0.20.1`

Commands:
```bash
mkdir -p trimmed/paired trimmed/unpaired
echo $SAMPLES | tr " " "\n" | parallel -j3 \
        fastp \
        --trim_poly_g \
        --average_qual 25 \
        --length_required 20 \
        --cut_front --cut_tail \
        --cut_mean_quality 25 --cut_window_size 10 \
        --trim_front1 5 --trim_front2 5 \
        --trim_tail1 5 --trim_tail2 5 \
        -i raw/{}_1.fastq.gz \
        -I raw/{}_2.fastq.gz \
        -o trimmed/paired/{}_1.fastq.gz \
        -O trimmed/paired/{}_2.fastq.gz \
        --unpaired1 trimmed/unpaired/{}_unpaired.fastq.gz
```

# Step 3: Assembly (MetaSPAdes)

Conda environment: `conda install -c bioconda spades=3.14.1`

Commands:
```bash

# Concatenate all the samples
cat trimmed/paired/*_1.fastq.gz > /tmp/all_R1.fastq.gz
cat trimmed/paired/*_2.fastq.gz > /tmp/all_R2.fastq.gz
cat trimmed/paired/*_unpaired.fastq.gz > /tmp/all_unpaired.fastq.gz

# Run metaspades
metaspades.py \
    -1 /tmp/all_R1.fastq.gz \
    -2 /tmp/all_R2.fastq.gz \
    -s /tmp/all_unpaired.fastq.gz \
    -o ms_results -t $THREADS -m $MEM
```

# Step 4: Coverage (bwa)

Conda environment: `conda install -c bioconda -c conda-forge parallel bwa=0.7.17`

Commands:
```bash
mkdir bwa && cd bwa

# Make bwa index
bwa index -p reference ../ms_results/scaffolds.fasta

#===========================#
#   File make_coverage.sh   #
#===========================#

reads_dir="../trimmed/paired"
sample="$1"

bwa mem -t $THREADS -a -M reference \
    $reads_dir/${sample}_1.fastq.gz \
    $reads_dir/${sample}_2.fastq.gz \
    | samtools view -@ $THREADS -bh \
    | samtools sort -@ $THREADS -o ${sample}.bam

samtools index -@ $THREADS ${sample}.bam ${sample}.bai

#===========================#
#      END OF FILE          #
#===========================#

# Run bwa mem on each sample
echo "015 117 250" | parallel -j $THREADS bash make_coverage.sh {}
```

# Step 5: Download the draft genomes from ONT assembly

Conda environment: `conda install entrez-direct` (also available with [docker](quay.io/biocontainers/entrez-direct))
BioProject: PRJNA529454

```bash
esearch -db nuccore -query "PRJNA529454" | efetch -format fasta > assembly_ont.fasta
```

# Step 6: De-duplicate ONT assembly

No conda environment --> [docker container](nakor/stampede-clustergenomes:0.9.0)
```bash
mkdir clustergenomes-ONT && cd clustergenomes-ONT
Cluster_genomes.pl -f assembly_ont.fasta -c 80 -i 95

# Outputs:
# - assembly_ont_95-80.clstr
# - assembly_ont_95-80.fna
# - assembly_ont-cover.csv
# - assembly_ont-nucmer.out.coords
# - assembly_ont-nucmer.out.delta
```

>> 1322/1880 contigs remaining

`conda install -c bioconda minimap2=2.17`

# Step 7: Align illumina contigs on ONT assembly
```
minimap2 clustergenomes-ONT/assembly_ont_95-80.fna ms_results/scaffolds.fasta -t $THREADS > alignments.tsv
```

# Step 8: Make contig assignment mapping
```python
import pandas as pd

min_id = 0.95
lin_len = 2048

cols = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", 
        "matches", "alnBlockLen", "quality", "unknown1", "unknown2", "unknown3", "unknown4", 
        "unknown5", "unknown6"]

aln = pd.read_csv("alignments.tsv", header=None, names=cols, sep='\t')
aln = aln[((aln.matches / aln.qlen) >= min_id) & (aln.qlen > min_len)]
best_matches = aln.groupby('qname').matches.idxmax()

best_matches['tname'] = aln.tname[aln.matches].tolist()
tname_mapping = pd.Series({ name: i for i,name in enumerate(best_matches['tname']) })
best_matches['tid'] = tname_mapping[best_matches['tname']].tolist()


```
