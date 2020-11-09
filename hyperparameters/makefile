INPUT ?= data/station-ALOHA
OUTPUT ?= results/$(notdir $(INPUT))

sim:
	mkdir -p $(OUTPUT) && \
	python clustering_opt.py run \
	    --threads 60 \
		--fasta $(INPUT)/assembly.fasta \
		--h5 $(INPUT)/coverage_contigs.h5 \
		--output $(OUTPUT) \
		--n-train 3000000 \
		--debug
