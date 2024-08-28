# 10X Multiome Scripts

Scripts for the 10X multiome analysis (RNA + ATAC) in trisomy and disomy foetal liver samples, as related to the following manuscript:

**Marderstein, A.R. et al. (2024). Single-cell multi-omics map of human foetal blood in Down's Syndrome. _Nature_.**

- **SCAVENGE** (includes the primary script for running SCAVENGE)
- **SCENT** (includes any scripts relevant to the SCENT analysis)
- **MIRA** (Python notebooks detailing the MIRA analysis led by Jon Bezney)
- **etc** (includes misc scripts)
- **andrew** is a folder containing draft scripts that is currently being kept for housekeeping purposes

## Table of Contents
- [Introduction](#introduction)
  - [Data Availability](#data-availability)
  - [Other Information](#other-information)
- [Scripts](#scripts)
  - [MIRA](#mira)
  - [Motifs](#motifs)
  - [SCENT](#scent)
    - [Preprocessing for SCENT Input](#preprocessing-for-scent-input)
    - [Running SCENT - Bash Parallel Submission](#running-scent---bash-parallel-submission)
    - [Running SCENT - R Parallelization](#running-scent---r-parallelization)
    - [Running SCENT - R Function](#running-scent---r-function)
    - [Postprocessing of SCENT Output](#postprocessing-of-scent-output)
    - [Peak-gene Links - Chromosomal Distribution (Histogram)](#peak-gene-links---chromosomal-distribution-histogram)
    - [Peak-gene Links - Number of Discoveries (Barplot)](#peak-gene-links---number-of-discoveries-barplot)
    - [Accessibility-by-Trisomy Interaction - Number of Discoveries (Barplot)](#accessibility-by-trisomy-interaction---number-of-discoveries-barplot)
    - [Summary Table](#summary-table)
    - [Peak-gene Links + RBC GWAS Enrichments (Forest Plot)](#peak-gene-links--rbc-gwas-enrichments-forest-plot)
    - [Misc Script for Analyses](#misc-script-for-analyses)
  - [SCAVENGE](#scavenge)
    - [Running SCAVENGE](#running-scavenge)
    - [Lifting Over Fine-map SNPs from hg19 to hg38](#lifting-over-fine-map-snps-from-hg19-to-hg38)
    - [HSC Branch Analysis](#hsc-branch-analysis)
    - [HSC Branch Plots](#hsc-branch-plots)
  - [Somatic](#somatic)
    - [Hasaart et al - Liftover hg19 to hg38](#hasaart-et-al---liftover-hg19-to-hg38)
  - [Other](#other)
    - [QC-related](#qc-related)
    - [Multiome Celltype Freq - Differential Abundance](#multiome-celltype-freq---differential-abundance)
    - [Multiome Celltype Freq - Comparison with scRNA](#multiome-celltype-freq---comparison-with-scrna)
    - [GATA1 Accessibility](#gata1-accessibility)
    - [Plot Peak Accessibility + Gene Expression Tracks for TFR2, TSPAN32](#plot-peak-accessibility--gene-expression-tracks-for-tfr2-tspan32)
    - [Pseudobulk Expression of HSCs in Large scRNA-seq Data](#pseudobulk-expression-of-hscs-in-large-scrna-seq-data)
    - [Intersect Peaks with Fine-Mapped GWAS SNPs](#intersect-peaks-with-fine-mapped-gwas-snps)
    - [TFR2 Overexpression Experiments (Statistical Analysis)](#tfr2-overexpression-experiments-statistical-analysis)
  - [Anndata to Seurat](#anndata-to-seurat)

## Introduction

### Data Availability

Input data for scripts are based on the datasets that have been deposited on ArrayExpress.

The following data has been deposited on ArrayExpress: 
- **scRNA-seq FASTQ raw data and CellRanger count matrices** (accession number E-MTAB-13067)
- **10x Visium FASTQ raw data, SpaceRanger count matrices, run summary metrics, and spatiality outputs** (E-MTAB-13062)
- **Multiome snRNA-seq and snATAC-seq FASTQ raw data, CellRanger ARC count matrices, and ATAC fragment files** (E-MTAB-13070). 

### Other Information
You will need to install packages that are listed in the header of scripts prior to running them.

Can't find code relevant to the analysis that you are interested in? Please look here first:
- **[GitLab repository for spatial transcriptomics and other scRNA-seq analyses](https://gitlab.com/cvejic-group/downsyndrome/)**
- **[GitHub repository for 10X multiome analyses](https://github.com/drewmard/t21_multiome)**

Please reach me at andrew.marderstein@gmail.com if there are questions about the analysis.

## Scripts

### MIRA

This folder contains 4 Jupyter notebooks that contain all of the analysis related to the MIRA and CELLRANK trajectories.

### Motifs

This folder contains scripts relevant to TF motif analyses.

### SCENT

#### Preprocessing for SCENT Input
- **pre_scent.R**: Uses assemble_parallel_files.R to create the SCENT input

#### Running SCENT - Bash Parallel Submission
- **run_scent_wrapper.sh**: Wrapper for running SCENT (using **run_scent.R**, which in turn uses **SCENT.R**)

#### Running SCENT - R Parallelization
- **run_scent.R**

#### Running SCENT - R Function
- **SCENT.R**

#### Postprocessing of SCENT Output
- **post_scent.R**: Compiles all SCENT output into a single file

#### Peak-gene Links - Chromosomal Distribution (Histogram)
- **peak_gene_summarize_plot.R**

#### Peak-gene Links - Number of Discoveries (Barplot)
- **multiome_t21_v_h_barplot.R**

#### Accessibility-by-Trisomy Interaction - Number of Discoveries (Barplot)
- **interaction_barplot.R**

#### Summary Table
- **create_scent_supp_table.R**

#### Peak-gene Links + RBC GWAS Enrichments (Forest Plot)
- **peak_gene_enrichment_plot.R**

#### Misc Script for Analyses
- **post_int.R**

### SCAVENGE

#### Running SCAVENGE
- **SCAVENGE.R**

#### Lifting Over Fine-map SNPs from hg19 to hg38
- **finemap_liftover_scavenge.R**

#### HSC Branch Analysis
- **scavenge_statistical_analysis.R**

#### HSC Branch Plots
- **TRS_by_branch_plots.R**

### Somatic

#### Hasaart et al - Liftover hg19 to hg38
- **hasaart_liftover.sh**

### Other

#### QC-related
- **qc_metrics_supptab.R**, **qc_metrics.R**

#### Multiome Celltype Freq - Differential Abundance
- **multiome_prop.R**

#### Multiome Celltype Freq - Comparison with scRNA
- **multiome_scrna_comp.R**

#### GATA1 Accessibility
- **GATA1_geneactivity.R**

#### Plot Peak Accessibility + Gene Expression Tracks for TFR2, TSPAN32
- **plot_region_of_interest.R**

#### Pseudobulk Expression of HSCs in Large scRNA-seq Data
- **HSC_pb.R**

#### Intersect Peaks with Fine-Mapped GWAS SNPs
- **rbc_peaks_overlap.sh**

#### TFR2 Overexpression Experiments (Statistical Analysis)
- **TFR2_experiment.R**

### Anndata to Seurat
- **anndata_to_seurat.R**
