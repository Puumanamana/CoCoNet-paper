ARGS ?= 
PHASE ?= clustering
INPUT ?= data/station-ALOHA
OUTPUT ?= results/$(notdir $(INPUT))


sim:
	mkdir -p $(OUTPUT) && \
	python $(PHASE)_opt.py run $(ARGS) \
	    --threads 60 \
		--fasta $(INPUT)/assembly.fasta \
		--h5 $(INPUT)/coverage_contigs.h5 \
		--output $(OUTPUT) \
		--n-train 3000000 \
		--debug

aloha:
	mkdir -p $(OUTPUT) && \
	python $(PHASE)_opt.py run $(ARGS) \
	    --threads 60 \
		--fasta $(INPUT)/*.fasta \
		--bam $(INPUT)/coverage/*.bam \
		--output $(OUTPUT) \
		--tlen-range 200 500 \
		--n-train 3000000 \
		--debug
