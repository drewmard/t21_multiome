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



summarize_motif_counts_thresholds <- function(input_data,threshold){
temp <-  input_data %>%
    filter_motif_data(log10_p_thresh = threshold) 
out <- assay(temp,'negLog10Padj') %>%
    data.frame() %>%
    mutate( max = case_when(disomy.down > ts21.up ~ 'disomy.down',TRUE~'ts21.up') ) %>%
    dplyr::count(max) %>%
    mutate(thresh = threshold) 
out
}



###### LOAD DATA
bioparam <- BiocParallel::SerialParam()
pwms <- getMatrixSet(JASPAR2022,opts = list(matrixtype = "PWM",tax_group = "vertebrates"))
ABC_all_path <- '/oak/stanford/groups/smontgom/epadhi/ts21_motif_analysis/ABC_data/HSCsts21_EnhancerPredictionsAllPutative.tsv.gz'
all_peaks <- classify_peaks(adj_thresh = .1) 
expression_DE <- classify_expression() 
EPD_promoter_annotations_DE_genes <- annotate_EPD_promoters_expression(expression_DE) 
EPD_ATAC_overlap_expression_annotated <- annotate_ATAC_peaks_EPD_promoters(all_peaks,EPD_promoter_annotations_DE_genes) %>% 
    annotate_ABC(ABC_all_path)
temp_filtered <- EPD_ATAC_overlap_expression_annotated %>%
    data.frame() %>%
    group_by(seqnames,start) %>%
    arrange(desc(abs(logFC_ATAC))) %>%
    mutate(rank = row_number()) %>%
    filter(rank == 1) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig <-temp_filtered %>%
    filter(adj_sig_ATAC == 'sig' & promoter == TRUE) %>%
    adjust_length() %>%
    wrap_motif_enrichment('direction_ATAC',pwms)
AP1_motif_list <- rowData(ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig)$motif.name %>%
    data.frame() %>%
    filter(str_detect(.,'JUN|FOS|ATF')) %>%
    pull(.)

assayNames(ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig)
pdf('temp_heatmap.pdf',width = 12,height = 7)
plotMotifHeatmaps(res %>% filter_motif_data(log10_p_thresh =1.3),cluster= TRUE,which = c('negLog10Padj','log2enr'),maxSig = 10,,show_seqlogo = TRUE)
dev.off()



res <- pull_motifs_monalisa_output(ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig,AP1_motif_list)
#ABC_vs_EPD_motif_enr <- EPD_ATAC_overlap_expression_annotated %>%
    #filter(!(promoter == FALSE & HSCsall == FALSE))  %>%
    #adjust_length() %>%
    #wrap_motif_enrichment('promoter',pwms)

#pdf('ABC_vs_EPD_heatmap.pdf',width = 12,height = 35)
#plotMotifHeatmaps(filt,cluster = TRUE,maxEnr = 2,which  = c('log2enr','negLog10Padj'),maxSig = 10)
#dev.off()


####### BINNED MOTIF ENRICHMENT, WITH BACKGROUND
# downsample 0 bin so that it is balanced with all 
# other bins and rejoin data 
zero_bin <- all_peaks %>%
    filter(bin == 0) %>%
    dplyr::sample_n(5000)
balanced_data <- binned_data %>%
    data.frame() %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE)

# gets sequences of the balanced data and runs monalisa across bins 
balanced_peask <- getSeq(BSgenome.Hsapiens.UCSC.hg38,balanced_data)
direction_binned <- calcBinnedMotifEnrR(seqs = balanced_peask,bins =as.factor(balanced_data$bin) ,pwmL =pwms)
saveRDS(direction_binned,'output_data/direction_binned_motif_enrich.RDS')

direction_binned_filtered <- direction_binned %>%
    filter_motif_data(log10_p_thresh=4) %>%
    create_summary_data

summarized_motif_data <- direction_binned_filtered %>%
    filter(!is.na(target)) %>% 
    rowwise() %>%
    mutate(max = max(X0,X1,X2)) %>%
    group_by(target) %>% 
    filter(max(max) == max) %>%
    ungroup() %>%
    dplyr::select(-max)

GATA_list <-summarized_motif_data %>%
    filter(str_detect(name,'GATA')) %>%
    pull(name)
chr21_list <- summarized_motif_data  %>%
    filter(chromosome_name == '21') %>%
    pull(name)


pdf('plots/direcion_binned_motif_enrich.pdf')
plotMotifHeatmaps(x = direction_binned[c(chr21_list,GATA_list),], which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE)
dev.off()



######## TS21 VS DISOMY ENRICHMENT, NO BACKGROUND
# filter to peaks that are either signifcantly different 
# in either ts21 or disomy and excludes nonsignifcant  sequences
ts21_vs_disomy_enr <-all_peaks %>%
    filter(adj_sig_ATAC == 'sig') %>%
    wrap_motif_enrichment('direction_ATAC',pwms)


# rerun motif analysis but now correcting for length
ts21_vs_disomy_enr_length_adjust <-all_peaks %>%
    filter(adj_sig_ATAC == 'sig') %>%
    adjust_length() %>%
    wrap_motif_enrichment('direction_ATAC',pwms)

saveRDS(ts21_vs_disomy_enr_length_adjust,'output_data/ts21_vs_disomy_enr_length_adjust.RDS')
ts21_vs_disomy_enr_length_adjust <- readRDS('../output_data/ts21_vs_disomy_enr_length_adjust.RDS')

saveRDS(ts21_vs_disomy_enr,'output_data/ts21_vs_disomy_no_background/ts21_vs_disomy_enr.RDS')
ts21_vs_disomy_enr <- readRDS('output_data/ts21_vs_disomy_no_background/ts21_vs_disomy_enr.RDS')

adjustment_comparison_pvals <- plot_adjustment_comparison(ts21_vs_disomy_enr,ts21_vs_disomy_enr_length_adjust,plot = FALSE)
adjustment_comparison_pvals <- plot_adjustment_comparison(ts21_vs_disomy_enr,ts21_vs_disomy_enr_length_adjust,plot = FALSE) + xlim(0,25) + ylim(0,25)
ggsave('plots/ts21_vs_disomy_adjustment_comparison.pdf')


motif_count_thresholds <- seq(1,20) %>%
    data.frame() %>%
    pmap_dfr(~summarize_motif_counts_thresholds(ts21_vs_disomy_enr_length_adjust,.)) %>%
    dplyr::rename('Group' = 'max') %>%
    ggplot(aes( x = thresh,y = n,group = Group,color = Group)) +
        geom_line() +
        geom_point() + 
        geom_hline(yintercept = 0,linetype = 'dashed',alpha = 0.6) +
        theme_classic() +
        ylab('Number of motifs') +
        xlab('-log10(p) Threshold') + 
        theme(legend.position = c('top'))
ggsave('plots/count_threshold.pdf')        



# gets list of motif IDs from the motif enrichment 
# analysis and subsets the JASPAR dataframe 
motif_list <- ts21_vs_disomy_enr_filtered %>%
    pull(motif)
ts21_PWMS <- getMatrixByID(JASPAR2022,ID = motif_list )
ts21_PWMS_motif_similarity <- motifSimilarity(ts21_PWMS) 


# plot heatmap for a subset of TFs
selected_list <- c('GATA5','GATA2','GATA1','GATA4','BACH1','RUNX1','NFE2','FOS::JUNB')
selected_TF_indexes <- data.frame(motif = rowData(ts21_vs_disomy_enr)$motif.name) %>%
    mutate(index = row_number()) %>%
    filter(toupper(motif) %in% selected_list) %>%
    pull(index)
selected_TFs <- ts21_vs_disomy_enr[selected_TF_indexes,]

pdf('plots/selected_TF_ts21_vs_disomy_enr.pdf',height = 7)
plot_dat <- ts21_vs_disomy_enr_length_adjust %>%
    pull_motifs_monalisa_output(selected_list)
plotMotifHeatmaps(x = selected_TFs, which.plots = c("log2enr"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE)
dev.off()




####PLOT THE CORRELATION MATRIX BETWEEN MOTIFS 
# gets features to plot on heatmap
maxwidth <- max(sapply(TFBSTools::Matrix(ts21_PWMS), ncol))
seqlogoGrobs <- lapply(ts21_PWMS, seqLogoGrob, xmax = maxwidth)
hmSeqlogo <- rowAnnotation(logo = annoSeqlogo(seqlogoGrobs, which = "row"),
                           annotation_width = unit(1.5, "inch"),
                            annotation_height = unit(1.5, "inch"),
                           show_annotation_name = FALSE)

# plots similarity matrix of ts21 results
pdf('ts21_disomy_motif_similarity.pdf',height = 12,width = 12)
Heatmap(ts21_PWMS_motif_similarity, 
        show_row_names = TRUE, row_names_gp = gpar(fontsize = 8),
        show_column_names = FALSE, column_names_gp = gpar(fontsize = 8),
        name = "Similarity", column_title = "Selected TFs",
        col = colorRamp2(c(0, 1), c("white", "red"))) 
        #right_annotation = hmSeqlogo)
dev.off()


#### COMPARE MOTIF ENRICHMENT TO DIFFERENTIAL EXPRESSION 
ts21_disomy_motif_enrich_all_results <- ts21_vs_disomy_enr_length_adjust %>%
    filter_motif_data(log10_p_thresh = 1.3) %>%
    create_summary_data(type = 'log2enr') 

motif_enrich_expression <- ts21_disomy_motif_enrich_all_results %>%
    mutate(name = toupper(name)) %>%
    left_join(expression_DE,by = c('name' = 'names')) %>%
    filter(!if_any(everything(),~is.na(.))) 

enrichment_correlation_plot <- motif_enrich_expression %>% 
    ggplot(aes(x = ts21.up,y = logFC_expr))  +
        geom_point(size = 2) + 
        xlab('Log2 Fold Change Expression') + 
        ylab('Motif enrichment') + 
        geom_smooth(method = 'lm') +
        theme_classic() +
        annotate('text',x =-1.2,y = 2.5,label = 'rho = 0.266 \np = 0.002') + 
        theme(axis.text = element_text(size = 15),axis.title = element_text(size = 17))
ggsave('motif_enrich_expression_LFC_correlation.pdf')

####PLOT THE NUMBER OF ENRICHED MOTIFS PER CHROMOSOME
ts21_vs_disomy_summarized <- ts21_vs_disomy_enr_length_adjust %>%
    filter_motif_data(log10_p_thresh = 1.3) %>%
    create_summary_data %>%
    data.frame() %>%
    filter(!is.na(target)) %>% 
    rowwise() %>%
    mutate(max = max(disomy.down,ts21.up)) %>%
    group_by(target) %>% 
    filter(max(max) == max) %>%
    ungroup() %>%
    dplyr::select(-max) 
motif_count_ts21_vs_disomy <- ts21_vs_disomy_summarized %>% 
    pivot_longer(!c(motif,input:name,chromosome_name),names_to = 'cluster')  %>%
    mutate(chromosome_name = case_when(chromosome_name == '21' ~ '21',TRUE ~'Other')) %>% 
    filter(value > 2) %>% 
    dplyr::count(cluster,chromosome_name) %>% 
    ggplot(aes(x = cluster,y = n,fill = chromosome_name)) +
        geom_col() +
        scale_fill_manual(values = c('goldenrod3','#A3A7AB')) +
        ylab('Motif Count') + 
        theme_classic()
ggsave('plots/motif_counts_ts21_vs_disomy.pdf',plot = motif_count_ts21_vs_disomy)


###### MOTIF ENRICHMENT OF PROMOTERS VS NON PROMOTERS 
# calculate baseline difference in motifs between promoters and non promoters 
# This is a bit confoudned by GC content and length so further down well 
promoter_motif_enr <- EPD_ATAC_overlap_expression_annotated %>%  data.frame() %>% 
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE) %>%
    wrap_motif_enrichment('promoter',pwms)
saveRDS(promoter_motif_enr,'output_data/promoter_analysis/EPD_ATAC_vs_ATAC_motif_enr.RDS')
promoter_motif_enr <- readRDS('../output_data/promoter_analysis/EPD_ATAC_vs_ATAC_motif_enr.RDS')

# plot GC content and length for promoters and non promoters  
pdf('plots/promoter_analysis/bin_diagnostics.pdf')
input_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38,EPD_ATAC_overlap_expression_annotated)
input_bins <- as.factor(EPD_ATAC_overlap_expression_annotated$promoter)
plotBinDiagnostics(input_seqs,input_bins,aspect = c("length"))
plotBinDiagnostics(input_seqs,input_bins,aspect = c("GCfrac"))
dev.off()


# adjust for lenth, rejection sampling doesnt work likely because the difference in the means between 
# promoters and non promoter fragments is too large so instead i resize fragments to have a uniform length 
length_adj_EPD_ATAC_overlap <- trim(resize(EPD_ATAC_overlap_expression_annotated, width = median(width(EPD_ATAC_overlap_expression_annotated)), fix = "center")) %>%
    GetGC(BSgenome.Hsapiens.UCSC.hg38)

# use matchRanges to create a beter matched range set for GC content for the non promoter sequences
matched_ranges <- nullranges::matchRanges(focal = length_adj_EPD_ATAC_overlap %>% filter(promoter == TRUE),
                              pool = length_adj_EPD_ATAC_overlap %>% filter(promoter == FALSE),
                              covar = ~ GC_content,method = 'nearest' ,replace = TRUE) 
# create granges object with promoter entries and the 
# generated matched ranges
matched_promoter_ranges <- length_adj_EPD_ATAC_overlap %>%
    filter(promoter == TRUE) %>%
    data.frame() %>%
    bind_rows(matched_ranges %>% data.frame()) %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE)

# run motif analysis on new ranges data that has adjusted the enhancers 
# for GC content and all sequences for length
promoter_motif_enr_matched_ranges <-  matched_promoter_ranges  %>%
    wrap_motif_enrichment('promoter',pwms)
saveRDS(promoter_motif_enr_matched_ranges,'output_data/promoter_analysis/EPD_ATAC_motif_enrich_matched_ranges.RDS')
promoter_motif_enr_matched_ranges <- readRDS('../output_data/promoter_analysis/EPD_ATAC_motif_enrich_matched_ranges.RDS')




# plot adjusted motif results for promoter analysis
pdf('plots/EPD_ATAC_promoter_vs_non_promoter_adjusted_motifs_temp.pdf',height = 25)
plot_dat<- promoter_motif_enr_matched_ranges %>%
    filter_motif_data(log10_p_thresh=8) 
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE,show_seqlogo = TRUE)
dev.off()

# plot adjusted motif results for promoter analysis
pdf('plots/EPD_ATAC_promoter_vs_non_promoter_adjusted_motifs_select_promoter_motifs.pdf',height = 15,width = 10)
plot_dat<-promoter_motif_enr_matched_ranges %>%
    pull_motifs_monalisa_output(promoter_motifs)
    #filter_motif_data(log10_p_thresh=8) 
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10)
dev.off()


pdf('plots/EPD_ATAC_promoter_vs_non_promoter_adjusted_motifs_select_enhancer_motifs.pdf',height = 15,width = 10)
plot_dat<-promoter_motif_enr_matched_ranges %>%
    pull_motifs_monalisa_output(enhancer_motifs)
    #filter_motif_data(log10_p_thresh=8) 
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10)
dev.off()




###### MOTIF ENRICHMENT OF DA PROMOTERS 
# runs motif enrichment on EPD promoters that are nominally DA and  
# then bins by direction and corrects for length 
ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig <-EPD_ATAC_overlap_expression_annotated %>%
    filter(nominal_sig_ATAC == 'sig' & promoter == TRUE) %>%
    adjust_length() %>%
    wrap_motif_enrichment('direction_ATAC',pwms)
saveRDS(ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig,'output_data/ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig.RDS')
ts21_vs_disomy_EPDpromoters_ATAC_length_adjusted_ATAC_direction<- readRDS('output_data/ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig.RDS')

# runs motif enrichment on EPD promoters that are adjusted signficiant  DA and  
# then bins by direction and corrects for length 
ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig <- EPD_ATAC_overlap_expression_annotated %>%
    filter(adj_sig_ATAC == 'sig' & promoter == TRUE) %>%
    adjust_length() %>%
    wrap_motif_enrichment('direction_ATAC',pwms)
saveRDS(ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig,'output_data/ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig.RDS')


# samples 4000 non DA sequences that are not promoters to add 
# to background sequences of the promoter analysis 
negative_seqs <- EPD_ATAC_overlap_expression_annotated %>%
    filter(promoter == FALSE & nominal_sig_ATAC == 'nonsig') %>%
    data.frame() %>%
    sample_n(4000) %>%
    mutate(direction_ATAC = 'nonsig')

# runs motif enrichment on EPD promoters that are adjusted signficiant  DA and  
# then bins by direction and corrects for length 
ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background <- EPD_ATAC_overlap_expression_annotated %>%
    filter(adj_sig_ATAC == 'sig' & promoter == TRUE) %>%
    data.frame() %>%
    bind_rows(negative_seqs) %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE) %>%
    adjust_length() %>%
    wrap_motif_enrichment('direction_ATAC',pwms)
saveRDS(ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background,'output_data/ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background.RDS')
ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background <- readRDS('output_data/ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background.RDS')
filtered_EPD_analysis <- ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background %>% filter_motif_data(log10_p_thresh = 1.3)
filtered_EPD_analysis %>% create_summary_data() %>% filter(chromosome_name == 21)

pdf('plots/ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background.pdf',height = 15)
plot_dat<-  pull_motifs_monalisa_output(ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background,filtered_AP1_JASPAR_IDs) %>%
    filter_motif_data(log10_p_thresh=1,type = 'negLog10Padj') 
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE,show_seqlogo = TRUE)
dev.off()


# Make heatmap of promoters binned by direction of differential expression
pdf('plots/ts21_vs_disomy_EPDpromoters_ATAC_direction_adj.pdf',height = 15)
plot_dat<-ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig %>%
    filter_motif_data(log10_p_thresh=1,type = 'negLog10Padj') 
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE)
dev.off()

AP1_motif_list <- rowData(ts21_vs_disomy_EPDpromoters_ATAC_length_adjusted_ATAC_direction_adj_sig)$motif.name %>%
    data.frame() %>%
    filter(str_detect(.,'JUN|FOS|ATF')) %>%
    pull(.)



# Make heatmap map of promoters binned direction of differential accessibility just for AP1 Motifs 
pdf('plots/EPD_ATAC_overlap_ATAC_binned_motif_res_AP1.pdf',height = 15)
plot_dat<-pull_motifs_monalisa_output(ts21_vs_disomy_EPDpromoters_ATAC_length_adjusted_ATAC_direction_adj_sig,filtered_AP1_JASPAR_IDs)
    #filter_motif_data(log10_p_thresh=1.3) 
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = all_ap1_clust, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE,show_seqlogo = TRUE)
dev.off()

selected_list <- c('GATA5','GATA2','GATA1','GATA4','BACH1','RUNX1','NFE2','FOS::JUNB')
pdf('plots/EPD_ATAC_overlap_ATAC_binned_motif_res_selected_TFs.pdf',height = 15)
plot_dat<-pull_motifs_monalisa_output(ts21_vs_disomy_EPDpromoters_ATAC_length_adjusted_ATAC_direction_adj_sig,selected_list) 
    #filter_motif_data(log10_p_thresh=1.3) 
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE,show_seqlogo = TRUE)
dev.off()


#######  ENHANCER MOTIF ANALYSIS - non ABC analysis 
### This section performs motitf analysis with promoters vs evberything not annotated as a promoter

# filters out enhancers and runs length adjusted binned motif enrichment 
enhancer_ts21_vs_disomy_motif_length_adjust <-EPD_ATAC_overlap_expression_annotated %>%
    filter(adj_sig_ATAC == 'sig' & promoter == FALSE) %>%
    adjust_length() %>%
    wrap_motif_enrichment('direction_ATAC',pwms)
saveRDS(enhancer_ts21_vs_disomy_motif_length_adjust,'output_data/enhancer_ts21_vs_disomy_motif_length_adjust.RDS')
enhancer_ts21_vs_disomy_motif_length_adjust <- readRDS('output_data/enhancer_ts21_vs_disomy_motif_length_adjust.RDS')
# plot adjusted motif results for enhancer analysis 

pdf('plots/enhancer_analysis/ATAC_enhancers_ts21_vs_disomy.pdf',height = 25)
plot_dat<-enhancer_ts21_vs_disomy_motif_length_adjust %>%
    filter_motif_data(log10_p_thresh=8) 
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE)
dev.off()

# plot adjusted motif results for enhancer analysis 
pdf('plots/ATAC_enhancers_ts21_vs_disomy_AP1.pdf',height = 12)
plot_dat<- pull_motifs_monalisa_output(enhancer_ts21_vs_disomy_motif_length_adjust,filtered_AP1_JASPAR_IDs) 
    #filter_motif_data(log10_p_thresh = 1.3)
plotMotifHeatmaps(x = plot_dat, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE,show_seqlogo = TRUE,cluster =all_ap1_clust )
dev.off()
####### ENHANCER MOTIF ANALYSIS - ABC
### This section performs a similar analysis to the one above but now using 
### ABC to define what an enhancer is to then run motif analysis stratified 
### by direction of differential accessibility

# filters out enhancers and runs length adjusted binned motif enrichment 
ABC_ts21_vs_disomy_motif_length_adjust <-EPD_ATAC_overlap_expression_annotated %>%
    filter(adj_sig_ATAC == 'sig' & HSCsts21 == TRUE & promoter == FALSE) %>%
    adjust_length() %>%
    wrap_motif_enrichment('direction_ATAC',pwms)
saveRDS(ABC_ts21_vs_disomy_motif_length_adjust,'output_data/ABC_ts21_vs_disomy_motif_length_adjust.RDS')
ABC_ts21_vs_disomy_motif_length_adjust <- readRDS('output_data/ABC_ts21_vs_disomy_motif_length_adjust.RDS')
filtered_ABC_res <-ABC_ts21_vs_disomy_motif_length_adjust %>% filter_motif_data(log10_p_thresh = 1.3)


filtered_ABC_res_summary_data <- create_summary_data(filtered_ABC_res)
TF_complexes <- filtered_ABC_res_summary_data %>%
    filter(is.na(input))
motif_name_id <- rowData(filtered_ABC_res)[,c('motif.id','motif.name')]
motif_name_id %>%
    data.frame() %>%
    filter(motif.id %in% TF_complexes$motif)


pdf('plots/ABC_ts21_vs_disomy_motif_direction_heatmap.pdf',height = 15)
plotMotifHeatmaps(filtered_ABC_res,which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE, show_seqlogo = TRUE)
dev.off()



####### PEAK GENE MOTIF ENRICH 
# READ IN PEAKS FROM ENHANCER TO GENE MAPPING FILE
peak_gene_df <- fread('/oak/stanford/groups/smontgom/epadhi/ts21_motif_analysis/DownSyndrome_HSC_PeakGeneSets/peak_gene.txt') %>%
    dplyr::select(peak,everything()) %>%
    separate(peak,sep = '-',into = c('chr','start','end'))

# create bins based on FDR of t21 term and 
# run motif enrichment and plot results  
ts21_beta_stratified_fdr <- peak_gene_df %>%
    #filter(!is.na(fdr_t21)) %>%
    mutate(bin = case_when(fdr_t21 < .1 ~1,TRUE ~ 2 )) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
ts21_beta_stratified_fdr_enrich <- wrap_motif_enrichment(ts21_beta_stratified_fdr,ts21_beta_stratified_fdr$bin,pwms)
ts21_beta_stratified_fdr_enrich_filtered <- ts21_beta_stratified_fdr_enrich %>%
        filter_motif_data(log10_p_thresh =2)
saveRDS(ts21_beta_stratified_fdr_enr,'/home/epadhi/epadhi/ts21_motif_analysis/analysis/output_data/E2G_analysis/ts21_beta_stratified_fdr.RDS')



pdf('plots/E2G_motif_analysis/e2G_fdr_ts21_stratified_enrichment.pdf',height = 15)
plotMotifHeatmaps(x = summar, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE)
dev.off()

# stratify peaks by pvalue and run motif enrichment 
# and plots results
ts21_beta_stratified_pval <- peak_gene_df %>%
    mutate(bin = case_when(boot_basic_p_t21 < .05 ~1,TRUE ~ 2 )) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
ts21_beta_stratified_pval_enr <- wrap_motif_enrichment(ts21_beta_stratified_pval,ts21_beta_stratified_pval$bin,pwms)
ts21_beta_stratified_pval_enrich_filtered <- ts21_beta_stratified_pval_enr %>%
        filter_motif_data(log10_p_thresh =2)
saveRDS(ts21_beta_stratified_pval_enr,'/home/epadhi/epadhi/ts21_motif_analysis/analysis/output_data/E2G_analysis/ts21_beta_stratified_pval.RDS')


pdf('plots/E2G_motif_analysis/e2G_pval_ts21_stratified_enrichment.pdf',height = 15)
plotMotifHeatmaps(x = ts21_beta_stratified_pval_enrich_filtered, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE)
dev.off()

# stratify results by pvalue on the disomy term from 
# the E2G model and run motif enrichment
disomy_beta_stratified_pval <- peak_gene_df %>%
    mutate(bin = case_when(boot_basic_p_H < .05 ~1,TRUE ~ 2 )) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
disomy_beta_stratified_pval_enr <- wrap_motif_enrichment(disomy_beta_stratified_pval,disomy_beta_stratified_pval$bin,pwms)
saveRDS(disomy_beta_stratified_pval_enr,'/home/epadhi/epadhi/ts21_motif_analysis/analysis/output_data/E2G_analysis/disomy_beta_stratified_pval.RDS')


disomy_beta_stratified_pval_enrich_filtered <- disomy_beta_stratified_pval_enr %>%
        filter_motif_data(log10_p_thresh =2)
pdf('plots/E2G_motif_analysis/e2G_pval_disomy_stratified_enrichment.pdf',height = 15)
plotMotifHeatmaps(x = disomy_beta_stratified_pval_enrich_filtered, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, maxEnr = 2, maxSig = 10,show_dendrogram = TRUE)
dev.off()


