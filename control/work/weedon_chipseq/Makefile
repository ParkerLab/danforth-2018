ANALYSIS_PATH=$(DANFORTH_HOME)/work/weedon_chipseq
ANALYSIS_NAME="WEEDON"

pipeline: $(ANALYSIS_PATH)/pipeline

$(ANALYSIS_PATH)/pipeline: commands
	@python commands

run: $(ANALYSIS_PATH)/pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
