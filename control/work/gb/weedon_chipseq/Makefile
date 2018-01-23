ANALYSIS_PATH=$(DANFORTH_HOME)/work/gb/weedon_chipseq
ANALYSIS_NAME="WN"

.PHONY: libraries 

libraries:
	@cd setup && bash commands && cd .. && ln -sf setup/libraries .

pipeline: $(ANALYSIS_PATH)/pipeline

$(ANALYSIS_PATH)/pipeline: commands libraries
	@python commands

run: $(ANALYSIS_PATH)/pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
