CONTROL=$(DANFORTH_HOME)/control

.PHONY: data setup atacseq atacseq_normalization ataqv_session differential_gene_expression differential_peaks master_peaks masui rnaseq rnaseq_normalization roadmap_chipseq roadmap_normalization weedon_chipseq weedon_normalization  h3k4me1_barplot kegg_volcano_plot rnaseq_heatmap atacseq_and_rnaseq_volcano_plots

data:
	@cd $(DANFORTH_HOME) && cd data && make all && cd $(DANFORTH_HOME)

setup:
	@Rscript -e 'install.packages("src/tss.tar.gz")'


## ANALYSES
atacseq:
	@cd $(CONTROL)/work/atacseq && make run && cd $(DANFORTH_HOME)

atacseq_normalization:
	@cd $(CONTROL)/work/gb/atacseq && make run && cd $(DANFORTH_HOME)

ataqv_session:
	cd $(CONTROL)/work/ataqv_session && bash commands && cd $(DANFORTH_HOME)

differential_gene_expression:
	@cd $(CONTROL)/work/differential_gene_expression && make run && cd $(DANFORTH_HOME)

differential_peaks:
	@cd $(CONTROL)/work/differential_peaks && make run && cd $(DANFORTH_HOME)

go:
	@cd $(CONTROL)/work/go && make run && cd $(DANFORTH_HOME)

master_peaks:
	@cd $(CONTROL)/work/master_peaks && make run && cd $(DANFORTH_HOME)

masui:
	@cd $(CONTROL)/work/masui_blat && make run && cd $(DANFORTH_HOME)

rnaseq:
	@cd $(CONTROL)/work/rnaseq && make run && cd $(DANFORTH_HOME)

rnaseq_normalization:
	@cd $(CONTROL)/work/gb/rnaseq && make run && cd $(DANFORTH_HOME)

roadmap_chipseq:
	@cd $(CONTROL)/work/roadmap_chipseq && make run && cd $(DANFORTH_HOME)

roadmap_normalization:
	@cd $(CONTROL)/work/gb/roadmap_chipseq && make run && cd $(DANFORTH_HOME)

weedon_chipseq:
	@cd $(CONTROL)/work/weedon_chipseq && make run && cd $(DANFORTH_HOME)

weedon_normalization:
	@cd $(CONTROL)/work/gb/weedon_chipseq && make run && cd $(DANFORTH_HOME)

## FIGURES
h3k4me1_barplot:
	@cd $(CONTROL)/figures/h3k4me1_barplot && bash commands && cd -

kegg_volcano_plot:
	@cd $(CONTROL)/figures/kegg_volcano_plot && bash commands && cd -

rnaseq_heatmap:
	@cd $(CONTROL)/figures/rnaseq_heatmap && bash commands && cd -

atacseq_and_rnaseq_volcano_plots:
	@cd $(CONTROL)/figures/atacseq_and_rnaseq_volcano_plots && bash commands && cd -

transcription_off_insertion:
	@cd $(CONTROL)/figures/transcription_off_insertion/commands && bash commands && cd -

## SUPPLEMENTAL TABLE
fpkm_table:
	@cd $(CONTROL)/supplemental_materials/fpkm_values && bash commands && cd -
