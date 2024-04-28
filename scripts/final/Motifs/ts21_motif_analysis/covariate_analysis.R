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
library(ggrepel)

source('../utils/monalisa_functions.R')
source('../utils/data_processing.R')
source('/oak/stanford/groups/smontgom/epadhi/1KG_cell_type_MPRA/analysis/utils/motif_analysis_funcs.R')





all_peaks <- classify_peaks(adj_thresh = .1) 
expression_DE <- classify_expression() 
EPD_promoter_annotations_DE_genes <- annotate_EPD_promoters_expression(expression_DE) 
EPD_ATAC_overlap_expression_annotated <- annotate_ATAC_peaks_EPD_promoters(all_peaks,EPD_promoter_annotations_DE_genes) 


####### PLOT DISTRIBION OF PEAK SIZES FOR TS21 VS DISOMY 
length_hist <- all_peaks %>%
    data.frame() %>%
    filter(nominal_sig_ATAC == 'sig') %>%
    ggplot(aes(x = log10(width),fill = direction_ATAC)) +
        geom_histogram(position = 'identity',alpha = 0.7,show.legend = FALSE) +
        xlab('log10(bp)') + 
        theme_classic() + 
        theme(axis.text = element_text(size = 15),axis.title = element_text(size = 17),legend.text = element_text(size=15),legend.title = element_text(size=17))

GC_content_dist <- GetGC(all_peaks,BSgenome.Hsapiens.UCSC.hg38) %>%
    data.frame() %>%
    ggplot(aes(x = GC_content,fill = direction_ATAC)) + 
        geom_histogram(position = 'identity',alpha = 0.7) +
        xlab('GC content') + 
        theme_classic() + 
        theme(axis.text = element_text(size = 15),axis.title = element_text(size = 17),legend.text = element_text(size=15),legend.title = element_text(size=17))

outplot <- length_hist + GC_content_dist
ggsave('ATAC_fragment_dist_ts21_vs_disomy.pdf',plot = outplot,width = 15)

# run motif enrichment analysis on significant vs nonsignifcant promoters  
# plot GC content and length for promoters and non promoters  
EPD_ATAC_overlap_expression_annotated %>%
    data.frame() %>%
    filter(nominal_sig_expr == 'sig') %>%
    plot_bin_diagnostics('direction_expr','plots/promoter_analysis/DE_promoters_ts21_vs_disomy_bin_QC.pdf')



