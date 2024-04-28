library(tidyverse)
library(monaLisa)
library(data.table)
library(plyranges)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(BiocParallel)
library(ComplexHeatmap)
library(patchwork)
source('../utils/monalisa_functions.R')
source('../utils/data_processing.R')
source('/oak/stanford/groups/smontgom/epadhi/1KG_cell_type_MPRA/analysis/utils/motif_analysis_funcs.R')





# load in all peaks from differential accesibility analysis and 
# annotate wether they are signifcantly up in ts21 (bin 1) or down in 
# ts21 (bin 2) or not significant (bin 0)
all_peaks <- classify_peaks(adj_thresh = .1) 


# load in differential expression results of ts21 vs disomy liver HSC and 
# annotate genes based on their significance (adjusted p val and unadjusted) and 
# the direction of effect
expression_DE <- classify_expression() 


# load in EPD promoter annotation and annotate if they belong 
# to a gene that is signifcantly DE between disomy and trisomy 
# and annotate the direction. Also widens the size of the EPD promoter  annotation 
EPD_promoter_annotations_DE_genes <- annotate_EPD_promoters_expression(expression_DE) 

# overlap ATAC seq fragments with EPD promoters to annotate 
# which regions are functioning as promoters and creates  granges object
EPD_ATAC_overlap_expression_annotated <- annotate_ATAC_peaks_EPD_promoters(all_peaks,EPD_promoter_annotations_DE_genes) 



####### QUANTIFY OVERLAP OF EPD PROMOTERS AND ATAC SEQ
# number of genes that are expressed in HSCs 
num_genes_DE <- expression_DE %>% nrow

# number of promoters for genes in EPD, number can be bigger 
# then total nubmer of genes because multiple promoters per gene 
num_promotes_EPD <- EPD_promoter_annotations_DE_genes %>% data.frame() %>% nrow

# get number of genes with an EPD promoter
num_genes_EPD <- EPD_promoter_annotations_DE_genes %>% data.frame()  %>%
    dplyr::select(geneSymbol) %>% distinct() %>% nrow

# number of EPD promoters that are overlapping ATAC peaks 
num_EPD_ATAC_overlap <- EPD_ATAC_overlap_expression_annotated %>% data.frame() %>% filter(promoter == TRUE) %>% nrow

# gets number of distinct genes that with an EPD - ATAC overlap 
num_genes_EPD_ATAC_overlap <- EPD_ATAC_overlap_expression_annotated %>% data.frame() %>% filter(promoter == TRUE) %>%
    dplyr::select(geneSymbol) %>% distinct() %>% nrow


# percent of genes with an EPD promoter
percent_genes_EPD <- num_genes_EPD/num_genes_DE
# percent of promoters that overlap that overlap an ATAC peak 
percent_EPD_ATAC <- num_EPD_ATAC_overlap/num_promotes_EPD
# percent of genes with a promoter that overlaps an ATAC peak
percent_genes_EPD_ATAC <- num_genes_EPD_ATAC_overlap/num_genes_DE

# calculate odds ratio for seeing a promoter with  a DE gene 
# if the promoter is also differntially accessible using  the nominal thresholds
RNA_ATAC_enrichment_nominal<- EPD_ATAC_overlap_expression_annotated %>%
    data.frame() %>%
    xtabs(~nominal_sig_ATAC + nominal_sig_expr,data = . ) %>% 
    data.frame() %>%
    pivot_wider(names_from = nominal_sig_expr,values_from = Freq) %>%
    column_to_rownames('nominal_sig_ATAC') %>%
    do(broom::tidy(fisher.test(.)))  %>%
    mutate(type = 'Nominal')

# calculate odds ratio for seeing a promoter with  a DE gene 
# if the promoter is also differntially accessible using the adj thresholds
RNA_ATAC_enrichments <- EPD_ATAC_overlap_expression_annotated %>%
    data.frame() %>%
    xtabs(~adj_sig_ATAC + adj_sig_expr,data = . ) %>% 
    data.frame() %>%
    pivot_wider(names_from = adj_sig_expr,values_from = Freq) %>%
    column_to_rownames('adj_sig_ATAC') %>%
    do(broom::tidy(fisher.test(.))) %>%
    mutate(type = 'Adjusted') %>%
    bind_rows(RNA_ATAC_enrichment_nominal)

enrichment_plot <- RNA_ATAC_enrichments %>%
    ggplot(aes(x = type, y = estimate)) + 
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = conf.low,ymax = conf.high)) +
        ylab('Enrichment of differential accessibile \n promoters linked to DE genes (OR)') + 
        xlab('Threshold Used') + 
        geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
        theme_classic() + 
        theme(axis.text = element_text(size = 15),axis.title= element_text(size = 17))

summary_data <- data.frame(percent=c(percent_genes_EPD,percent_EPD_ATAC,percent_genes_EPD_ATAC),number = c(num_genes_EPD,num_EPD_ATAC_overlap,num_genes_EPD_ATAC_overlap),name = c('Genes with 1 or \n more EPD promoters','EPD promoters \n overlap with ATAC','Genes with an \n EPD ATAC overlap ')) 
overlap_summary <- summary_data %>%
    ggplot(aes(x = name,y = percent)) +
        geom_col() +
        geom_text(aes(label = number),vjust = -.3) +
        theme_classic() + 
        xlab('') + 
        ylab('Percent') +
        theme(axis.text = element_text(size = 15),axis.title= element_text(size = 17))

combined_promoter_overlap_plot <- overlap_summary+ enrichment_plot 
ggsave('plots/promoter_analysis/EPD_overlap_summary.pdf',plot = combined_promoter_overlap_plot,width = 15)


