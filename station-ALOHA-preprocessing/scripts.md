# Download data: accessions: 
- SRR5677468 (15m)
- SRR8811962 (117m)
- SRR8811963 (250m)

# Global variables
```bash
export SAMPLES="SRR5677468 SRR8811962 SRR8811963"
export THREADS=10
export MEM=500G
```

## Conda environment: `conda install -c bioconda sra-tools=2.10.8`
## Commands:
```bash
mkdir raw
fasterq-dump --threads $THREAD --mem $MEM --split-files $SAMPLES --outdir raw
for f in $(ls raw/*.fastq); do
    gzip $f
done
```

# Trimming with fastp
## Parameters
    - Crop reads at beginning and end: 5 bp
    - Trimming with sliding window of size 10 when quality < 25
    - Minimum read length after trimming: 20 bp
    - Trim poly-G sequences
## Conda environment: `conda install -c bioconda -c conda-forge parallel fastp=0.20.1`
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

# Assembly with MetaSPAdes
## Conda environment: `conda install -c bioconda spades=3.14.1`
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

# Coverage with bwa mem
## Conda environment: `conda install -c bioconda -c conda-forge parallel bwa=0.7.17`
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
