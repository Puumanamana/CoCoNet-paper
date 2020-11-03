DB ?= refseq-viral.fasta

SHELL := /bin/bash

main:
	wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{1,2}.1.genomic.fna.gz \
	&& zcat viral.*.gz > refseq-viral.fasta \
	&& rm -f viral*.gz \
	&& python clean_db.py --fasta $(DB) --output refseq-viral-ACGT.fasta \
	&& python plot_composition_separation.py --db refseq-viral-ACGT.fasta
