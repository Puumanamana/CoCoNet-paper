#!/usr/bin/env bash

wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{1,2}.1.genomic.fna.gz \
    && zcat viral*.gz > refseq-viral.fasta \
    && rm -f viral*.gz
