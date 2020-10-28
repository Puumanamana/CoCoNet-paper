FOLDER ?= data/camisim_2000-genomes_15-samples_10X_6
SUFFIX = $(notdir $(FOLDER))

sim:
	echo "======== Processing $(FOLDER) ========" && \
	python main.py run \
	    --fasta $(FOLDER)/assembly.fasta \
	    --h5 $(FOLDER)/coverage_contigs.h5 \
	    --output results/output-$(SUFFIX) \
	    --features coverage

SA:
	echo "======== Processing $(FOLDER) ========" && \
	python main.py run \
	    --fasta $(FOLDER)/*.fasta \
	    --bam $(FOLDER)/*.bam \
	    --output results/output-$(SUFFIX) \
	    --features coverage
