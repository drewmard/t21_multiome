# 10X multiome scripts

Scripts for the 10X multiome analysis (RNA + ATAC) in trisomy and disomy foetal liver samples, as related to the following manuscript:

Marderstein, A.R. et al. (2023). Single-cell multi-omics map of human foetal blood in Down's Syndrome. biorxiv.

- **SCAVENGE** (includes the primary script for running SCAVENGE)
- **SCENT** (includes any scripts relevant to the SCENT analysis)
- **MIRA** (Python notebooks detailing the MIRA analysis led by Jon Bezney)
- **etc** (includes misc scripts)
- **andrew** is a folder containing draft scripts that is currently being kept for housekeeping purposes

### Other info

Input data for scripts are based on the datasets that have been deposited on ArrayExpress.

Install packages that are listed in the header of scripts prior to running them.

Please reach out if there are questions about the analysis: amarder@stanford.edu.

## Scripts



### MIRA

This folder contains 4 jupyter notebooks that contain all of the analysis related to the MIRA and CELLRANK trajectories.




### motif




### SCENT

#### Preprocessing for SCENT input
pre_scent.R (Uses assemble_parallel_files.R to create the SCENT input)

#### Running SCENT - bash parallel submission
run_scent_wrapper.sh (wrapper for running SCENT (using **run_scent.R**, which in turn uses **SCENT.R**)

#### Running SCENT - R parallelization
run_scent.R

#### Running SCENT - R function
SCENT.R

#### Postprocessing of SCENT output (compiles all SCENT output into a single file)
post_scent.R

#### Peak-gene links - chromosomal distribution (histogram)
peak_gene_summarize_plot.R

#### Peak-gene links - number of discoveries (barplot)
multiome_t21_v_h_barplot.R

#### Accessibility-by-trisomy interaction - number of discoveries (barplot)
interaction_barplot.R

#### summary table
create_scent_supp_table.R

#### Peak-gene links + RBC GWAS enrichments (forest plot)
peak_gene_enrichment_plot.R

#### Misc script for analyses
post_int.R





### SCAVENGE

#### Running SCAVENGE
SCAVENGE.R

#### Lifting over fine-map SNPs from hg19 to hg38
finemap_liftover_scavenge.R

#### HSC branch analysis
scavenge_statistical_analysis.R

#### HSC branch plots
TRS_by_branch_plots.R



### Somatic

#### Hasaart et al - liftover hg19 to hg38
hasaart_liftover.sh




### Other

#### QC-related
qc_metrics_supptab.R, qc_metrics.R
#### multiome celltype freq - differential abundance
multiome_prop.R

#### multiome celltype freq - comparison with scRNA
multiome_scrna_comp.R

#### GATA1 accessibility
GATA1_geneactivity.R

#### Plot peak accessibility + gene expression tracks for TFR2, TSPAN32
plot_region_of_interest.R

#### Pseudobulk expression of HSCs in large scRNA-seq data
HSC_pb.R

#### Intersect peaks with fine-mapped GWAS SNPs
rbc_peaks_overlap.sh

#### TFR2 overexpression experiments (statistical analysis)
TFR2_experiment.R



### 

#### anndata --> seurat
anndata_to_seurat.R
