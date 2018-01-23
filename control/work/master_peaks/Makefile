ANALYSIS_PATH=${DANFORTH_HOME}/work/master_peaks
ANALYSIS_NAME="MP"

clean:
	@rm -rf $(ANALYSIS_PATH)

libraries:
	@cd setup && bash commands && cd .. && ln -sf setup/libraries .

.PHONY: libraries clean

pipeline: $(ANALYSIS_PATH)/pipeline

$(ANALYSIS_PATH)/pipeline: commands libraries
	@python commands

run: $(ANALYSIS_PATH)/pipeline
	@cd $(ANALYSIS_PATH) && drmr -j $(ANALYSIS_NAME) pipeline && cd -
