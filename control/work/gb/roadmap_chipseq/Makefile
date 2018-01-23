ANALYSIS_PATH=$(DANFORTH_HOME)/work/gb/roadmap_chipseq
ANALYSIS_NAME="gbC"

.PHONY: libraries

libraries:
	@cd setup && bash commands && cd .. && ln -sf setup/libraries .

pipeline: commands libraries
	@python commands

run: pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
