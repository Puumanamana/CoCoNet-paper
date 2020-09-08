opt:
	mkdir -p sherpa-results-$(FOLDER)
	python main.py \
		--fasta $(FOLDER)/assembly.fasta \
		--h5 $(FOLDER)/coverage_contigs_gt2048.h5 \
		--output test-$(FOLDER)

aloha:
	mkdir -p sherpa-results-$(FOLDER)
	python main.py \
		--fasta $(FOLDER)/assembly.fasta \
		--bam $(FOLDER)/*.bam \
		--output test-$(FOLDER) \
		--fl-range 200 500

aloha_old:
	mkdir -p sherpa-results-$(FOLDER)
	python old.py \
		--fasta $(FOLDER)/assembly.fasta \
		--coverage $(FOLDER)/*.bam \
		--output test-$(FOLDER)
