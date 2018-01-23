ANALYSIS_PATH=$(DANFORTH_HOME)/work/atacseq
ANALYSIS_NAME="ATAC"

libraries:
	@cd setup && bash commands && cd .. && ln -sf setup/libraries .

.PHONY: libraries

pipeline: $(ANALYSIS_PATH)/pipeline

$(ANALYSIS_PATH)/pipeline: commands libraries
	@python commands

run: $(ANALYSIS_PATH)/pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
