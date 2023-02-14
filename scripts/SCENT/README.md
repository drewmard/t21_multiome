# SCENT

## Important scripts to run SCENT:

#### pre_scent.R

Uses **assemble_parallel_files.R** to create the input for SCENT analysis.

#### run_scent_wrapper.sh

This is a wrapper for running the SCENT analysis (using **run_scent.R**, which in turn uses **SCENT.R**).

#### post_scent.R

This compiles all SCENT output into a single file.

### Downstream analysis:

**create_scent_supp_table.R**: create SCENT supplementary table

**peak_gene_enrichment_plot.R**: make peak-gene enrichment forest plot w/ RBC GWAS

**plot_region_of_interest.R**: plot peak accessibility + gene expression tracks

**peak_gene_summarize_plot.R**: chromosomal distribution of peak-gene links

**HSC_pb.R**: pseudobulk plots of genes

