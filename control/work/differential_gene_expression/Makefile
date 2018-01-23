ANALYSIS_PATH=$(DANFORTH_HOME)/work/differential_gene_expression
ANALYSIS_NAME="DGE"

.PHONY: libraries

libraries:
	@cd setup && bash commands && cd .. && ln -sf setup/libraries .

pipeline: libraries
	@python commands

run: pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
