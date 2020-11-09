FOLDER ?= results/station-ALOHA

BAM := $(wildcard results/station-ALOHA/illumina-coverage/*.bam)

dataset:
	nextflow run nf-station-aloha -with-singularity --run_assembly --outdir $(FOLDER)

tlen: $(BAM)
$(BAM):
	samtools view $(@) \
	| awk '{print sqrt($$9*$$9)}' \
	| sort -n | uniq -c \
	| awk '{print $$1","$$2}' \
	> results/tlen-$(shell basename $(@) .bam).csv
.PHONY: tlen $(BAM)

plot:
	./plot_template-length.py
