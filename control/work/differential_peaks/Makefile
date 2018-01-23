ANALYSIS_PATH=$(DANFORTH_HOME)/work/differential_peaks
ANALYSIS_NAME="DP"

.PHONY: libraries

libraries:
	@cd setup && bash commands && cd .. && ln -sf setup/libraries .

pipeline: libraries commands
	@python commands

run: pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
