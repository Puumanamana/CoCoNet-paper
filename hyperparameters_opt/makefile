FOLDER ?= 500_3_4_0
OUTDIR ?= test-$(FOLDER)
ARGS ?= 

sim:
	python main.py run $(ARGS) \
	    --threads 10 \
		--fasta $(FOLDER)/assembly.fasta \
		--h5 $(FOLDER)/coverage_contigs.h5 \
		--output $(OUTDIR)

aloha:
	python main.py run $(ARGS) \
	    --threads 10 \
		--fasta $(FOLDER)/assembly.fasta \
		--bam $(FOLDER)/*.bam \
		--output $(OUTDIR) \
		--tlen-range 200 500
