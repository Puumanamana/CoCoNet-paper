FOLDER ?= 2000_10_15_0

sim:
	echo "======== Processing $(FOLDER) ========" && \
	python main.py run \
	    --fasta input_data/$(FOLDER)/assembly.fasta \
	    --h5 input_data/$(FOLDER)/coverage_contigs.h5 \
	    --output workspace/output-$(FOLDER) \
	    --features coverage

SA:
	echo "======== Processing $(FOLDER) ========" && \
	python main.py run \
	    --fasta input_data/$(FOLDER)/*.fasta \
	    --bam input_data/$(FOLDER)/*.bam \
	    --output workspace/output-$(FOLDER) \
	    --features coverage
