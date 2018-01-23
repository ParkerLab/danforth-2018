# 2018 Danforth analysis

This repository contains scripts to download data, run analyses, and recreate many of the figures from our 2018 Danforth analysis.

## General organization
This has been tested on a Linux server. Pipelines are meant to be used with drmr (can be downloaded from the Parker Lab GitHub). The general organization of the repo is as follows:

1. `src`, `bin`, and `sample_information` directories should be self-explanatory.
2. `data` directory contains our raw data and other 'generic' data such as fasta files, bwa indices, etc. (some elements are distributed in the repository itself, some are created/downloaded in the steps outlined below).
3. `control` directory contains the actual analysis scripts (setting up subdirectories in the `work` directory, and in many cases producing `pipeline` files that end up in the `work` directory).
4. `sw` contains some third-party software (RNA-Enrich).
5. The `work` directory will be created as the `make` commands outlined below are run. Some pipeline files will be printed into there, and the actual analysis results will be put in there (e.g. bam files, DESeq2 results, etc).
6. The `figures` directory will be created as the `make` commands outlined below are run. If figures are recreated using the `make` commands outlined later, they will appear in here (along with some necessary intermediate files).

There is a `Makefile` in the top level of the repository as well.  All the `make` commands mentioned below refer to this `Makefile` and should therefore be run from the top directory.

## Dependencies not included in repo
As stated, many of the pipelines for this analysis utilize `drmr` to submit jobs to a resource manager (we use SLURM). Therefore, `drmr` will need to be downloaded from the Parker Lab GitHub. Also, python (we used v. 2.7.13) and R (v. 3.3.3) will need to be present on the system.

### Tools run from the command line
These tools must be in your $PATH:
1. cta (v. 0.1.2; can be downloaded from the Parker Lab GitHub)
2. fastqc (v. 0.11.5)
3. bwa (v. 0.7.15-r1140)
4. picard (v. 2.8.1)
5. ataqv (v. 1.0; can be downloaded from the Parker Lab GitHub)
6. STAR (v. 2.5.2b)
6. QoRTs (v. 1.0.7)
6. bnMapper.py
6. macs2 (v. 2.1.1.20160309)
6. samtools (v. 1.3.1, using htslib 1.3.2)
6. bedtools (v. 2.26.0)
6. blat (v. 36x2)
6. SRA toolkit (for fastq-dump; v. 2.8.1)
6. fastx_trimmer (FASTX Toolkit v. 0.0.14)
6. mysql (for querying UCSC; v. 15.1 Distrib 10.1.26-MariaDB, for debian-linux-gnu (x86_64) using readline 5.2)
6. phantompeakqualtools (for ChIP-seq QC; v. 2.0)

### R packages
1. DESeq2 (v. 1.14.1)
2. dplyr (v. 0.7.4)
3. tidyr (v. 0.7.0)
10. ggplot2 (v. 2.2.1)
10. ggrepel (v. 0.6.5)
10. optparse (v. 1.4.4)
10. vsn (v. 3.42.3)
10. biomaRt (v. 2.30.0)
10. cowplot (v. 0.8.0)
10. QoRTs (v. 1.1.8)

### Additional python modules
1. pybigwig (v. 0.3.9)
1. pysam (v. 0.11.2.1)
1. numpy (v. 1.13.1)

## Setup
1. Clone the repository
2. Set your environmental variable $DANFORTH_HOME to the repository path (of course, you may want to add each of these `export` commands to your `.bashrc` so that you don't need to re-set them each time you log out and back in to the server):
```bash
export DANFORTH_HOME='/path/to/repo'
```
2. Add the src directory to your $PYTHONPATH and the bin directory to your $PATH:
```bash
export PYTHONPATH="$PYTHONPATH:${DANFORTH_HOME}/src"
export PATH="$PATH:${DANFORTH_HOME}/bin"
```
2. If interested in the ChIP-seq analyses, set a variable representing the path to phantompeakqualtool's `run_spp.R` script:
```bash
export RUN_SPP_PATH="/path/to/run_spp.r"
```
2. If interested in the RNA-seq analyses, set variables representing the path to QoRTs files:
```bash
export QORTS_JAR="/path/to/QoRTs.jar"
export QORTS_GEN_MULTI_QC="/path/to/qortsGenMultiQC.R"
```
and set a variable to indicate which genome the RNA-seq reads should be mapped to (mm9 for all analyses discussed in the manuscript; this variable is present only because it is needed in the case that one wants to recreate Fig. S1):
```bash
export DANFORTH_RNASEQ_GENOME="mm9"
```
2. Install a small included R package for linking peaks to their nearest TSS:
```bash
make setup
```
3. Prepare the `data` directory, containing all generic data (e.g., fasta files, bwa indices, etc). When the GEO repository becomes public, this will also download our ATAC-seq/RNA-seq fasta files from GEO (for now, this will obviously not happen):
```bash
make data
```
This will run most commands necessary to set up the data/ directory immediately in a consecutive fashion; however, because bwa indices take several hours to put together, it will submit jobs to drmr for the index creation (job names BWAINDEX). You should wait for those jobs to finish before proceeding, with the exception of processing the RNA-seq data as that does not use BWA.

## Running analyses

Now the pipelines can be run. Order is somewhat important, as some pipelines depend on others. Such dependences are noted below.

### ATAC-seq processing and differential peak calling

To carry out the primary processing of the ATAC-seq data (adapter trimming, mapping, duplicate removal, filtering, peak calling, creating ataqv sessions), one can just run (from the top level of the repository):
```bash
make atacseq
```
The results of this analysis will end up in the `${DANFORTH_HOME}/work/atacseq` directory. Once this pipeline has finished running, one can merge the peak calls to get the set of 'master peaks', and determine the number of reads that fall within each of these master peaks for each ATAC-seq experiment (this information will be used in the differential peak calling):
```bash
make master_peaks  # requires atacseq
```
The output from that will be in the `${DANFORTH_HOME}/work/master_peaks` directory. After this has completed, one can perform the differential peak calling:

```bash
make differential_peaks  # requires master_peaks
```
The output will be in the `${DANFORTH_HOME}/work/differential_peaks` directory.

To generate an ataqv session for the ATAC-seq data, run:
```bash
make ataqv_session
```

The session will be created in the `${DANFORTH_HOME}/work/ataqv_session` directory.

To create the signal tracks (bigwig files) that can be used to visually compare samples, one must normalize the signal for sequencing depth. This can be done by running:

```bash
make atacseq_normalization  # requires atacseq
```
The output will be in the `${DANFORTH_HOME}/work/gb/atacseq` directory.

### RNA-seq processing, differential gene expression analysis, and KEGG pathway enrichment analysis:

To carry out the primary processing of the RNA-seq data, one must run:
```bash
make rnaseq  # primary RNA-seq processing pipeline
```

The output will be in the `${DANFORTH_HOME}/work/rnaseq` directory. Once this has completed, the differential peak calling can be run:

```bash
make differential_gene_expression  # requires rnaseq
```

The output will be in the `${DANFORTH_HOME}/work/differential_gene_expression` directory. After this has finished, one can run the KEGG pathway enrichment analysis:
```bash
make go  # requires differential gene expression
```
The output will be in the `${DANFORTH_HOME}/work/go` directory.

To create the signal tracks (bigwig files) that can be used to visually compare samples, one must normalize the signal for sequencing depth (also, separate the signal by strand). This can be done by running:

```bash
make rnaseq_normalization  # requires rnaseq
```
The output will be in the `${DANFORTH_HOME}/work/gb/rnaseq` directory.

### Download and processing of the Weedon et al ChIP-seq data
Note that this will download and process more than just the H3K4me1 ChIP-seq data, but we only utilize the H3K4me1 data.

To do the downloading and processing:
```bash
make weedon_chipseq
```

The output will be in the `${DANFORTH_HOME}/work/weedon_chipseq` directory. In order to allow for comparisons with other experiments, the signal needs to be normalized for sequencing depth. This can be done by:
```bash
make weedon_normalization  # normalizes the bigwigs created by weedon_chipseq
```

The output will be in the `${DANFORTH_HOME}/work/gb/weedon_chipseq` directory.

### Download and processing of the Roadmap Epigenomics data

To download and process the Roadmap Epigenomics ChIP-seq data, one must run:
```bash
make roadmap_chipseq
```
The output will be in the `${DANFORTH_HOME}/work/roadmap_chipseq` directory. To normalize the signal from these experiments in order to allow for cross-tissue comparisons, once the above pipeline is complete you must run:

```bash
make roadmap_normalization
```
The output will be in the `${DANFORTH_HOME}/work/gb/roadmap_chipseq` directory.

### Determining locations of the Ptf1a binding sites from Masui et al:
```bash
make masui
```

## Re-creating figures

### RNA-seq volcano plot and ATAC-seq volcano plot
To create this figure, one must first have recreated the ATAC-seq data processing, master peak creation, and differential peak analysis outlined above, as well as the RNA-seq data processing and differential gene expression calling. Once these pipelines have been successfully run, one can create these volcano plots using the command:
```bash
make atacseq_and_rnaseq_volcano_plots
```

### H3K4me1 barplot
To create this figure, the Weedon et al data must have been processed and normalized, and the Roadmap Epigenomics data must have been processed and normalized. Also, because the signal refers to the signal over the region orthologous to the differential Gm13344 peak, the ATAC-seq data must have been processed and differential peak calling performed. Once this has been done, this figure can be created by running:

```bash
make h3k4me1_barplot
```

### RNA-seq heatmap
To create this figure, one must have processed the RNA-seq data and run the differential gene expression analysis. Once these pipelines have successfully completed, you can re-create this figures with the command:

```bash
make rnaseq_heatmap
```

### KEGG enrichment volcano plot
To create this figure, the KEGG enrichment analysis needs to have been run. Once this is done, the KEGG plot can be created using the command:
```bash
make kegg_volcano_plot
```

### Plot of RNA-seq signal at/near insertion (Fig S1)
To create this figure, the RNA-seq processing (`make rnaseq`) should be run with the `$DANFORTH_RNASEQ_GENOME` environmental variable set to `"danforth"` rather than `"mm9"`. Once this is done, run `make rnaseq_normalization` (keeping `$DANFORTH_RNASEQ_GENOME` set to `"danforth"`), and then run:
```bash
make transcription_off_insertion
```

## Re-creating table of FPKM values (in supplement):
To generate this table, run:
```bash
make fpkm_table
```

