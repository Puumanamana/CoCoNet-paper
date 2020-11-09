FOLDER ?= data/camisim_2000-genomes_15-samples_10X_6
SUFFIX = $(notdir $(FOLDER))

sim:
	echo "======== Processing $(FOLDER) ========" && \
	python main.py run --threads 50 \
	    --fasta $(FOLDER)/assembly.fasta \
	    --h5 $(FOLDER)/coverage_contigs.h5 \
	    --output results/output-$(SUFFIX) \
		--features coverage
plot:
	python plot.py results/scores.csv
