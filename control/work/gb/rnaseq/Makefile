ANALYSIS_PATH=$(DANFORTH_HOME)/work/gb/rnaseq
ANALYSIS_NAME="gbR"

.PHONY: libraries

libraries:
	@cd setup && bash commands && cd .. && ln -sf setup/libraries .

pipeline: $(ANALYSIS_PATH)/pipeline

$(ANALYSIS_PATH)/pipeline: commands libraries
	@python commands

run: pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
