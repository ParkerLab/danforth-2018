CONTROL_DIR=${DANFORTH_HOME}/control/data

# TODO:sequencing should be added in, once the data can be downloaded from GEO
all: chrom_sizes fasta bwa mappability tss gtf entrez_to_ensembl chain


.PHONY: clean chain chrom_sizes insertion mappability transcript_lengths tss sequencing entrez_to_ensembl

clean:
	@rm -rf ${DANFORTH_HOME}/data/bwa
	@rm -rf ${DANFORTH_HOME}/data/chain
	@rm -rf ${DANFORTH_HOME}/data/chrom_sizes
	@rm -rf ${DANFORTH_HOME}/data/entrez_to_ensembl
	@rm -rf ${DANFORTH_HOME}/data/fasta
	@rm -rf ${DANFORTH_HOME}/data/gtf
	@rm -rf ${DANFORTH_HOME}/data/mappability
	@rm -rf ${DANFORTH_HOME}/data/tss

bwa: fasta
	@python $(CONTROL_DIR)/bwa.commands
	@drmr -j BWAINDEX $(DANFORTH_HOME)/data/bwa/pipeline

chain:
	@bash $(CONTROL_DIR)/chain.commands

chrom_sizes:
	@bash $(CONTROL_DIR)/chrom_sizes.commands

fasta: insertion
	@bash $(CONTROL_DIR)/fasta.commands

mappability:
	@bash $(CONTROL_DIR)/mappability.commands

transcript_lengths:
	@bash $(CONTROL_DIR)/transcript_lengths.commands

tss:
	@bash $(CONTROL_DIR)/tss.commands

sequencing:
	@bash $(CONTROL_DIR)/sequencing.commands

gtf:
	@bash $(CONTROL_DIR)/gtf.commands

entrez_to_ensembl:
	@bash $(CONTROL_DIR)/entrez_to_ensembl.commands
