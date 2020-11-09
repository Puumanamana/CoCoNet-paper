N := 100 500 1000 5000 10000 50000 100000
R ?= 5

sim: $(N)
$(N):
	python simulation.py --db refseq-viral-ACGT.fasta --n-samples 5 --min-genome-len 3000 --n-contigs $(@) --replicate $(R)
.PHONY: sim $(N)
