ANALYSIS_PATH=$(DANFORTH_HOME)/work/roadmap_chipseq
ANALYSIS_NAME="RCHIP"

.PHONY: libraries

libraries:
	@cd setup && bash commands && cd .. && ln -sf setup/libraries libraries

pipeline: $(ANALYSIS_PATH)/pipeline

$(ANALYSIS_PATH)/pipeline: commands libraries
	@python commands

run: pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
