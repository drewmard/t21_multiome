get_motif_correlations <- function(input_motif_list,plot = TRUE){
motif_dat <- getMatrixByID(JASPAR2022,ID = input_motif_list)
motif_correlations <-  motifSimilarity(motif_dat) %>%
    data.matrix()
if (plot == TRUE){
maxwidth <- max(sapply(TFBSTools::Matrix(motif_dat), ncol))
seqlogoGrobs <- lapply(motif_dat, seqLogoGrob,xmax = maxwidth)
hmSeqlogo <- rowAnnotation(logo = annoSeqlogo(seqlogoGrobs, which = "row"),
                           annotation_width = unit(2, "inch"), 
                           show_annotation_name = FALSE
)
outplot <- Heatmap(motif_correlations,right_annotation = hmSeqlogo,show_column_names = FALSE)
outplot
}
else{
rownames(motif_correlations) <- input_motif_list 
motif_correlations
}
}




wrap_motif_enrichment <- function(input_df,bins,pwms){
require(monaLisa)
require(BSgenome.Hsapiens.UCSC.hg38)
input_bins <- input_df %>%
        data.frame() %>%
        pull({bins}) 
print(paste0('Number of input sequences: ', length(input_bins)))
peak_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38,input_df)
monalisa_res <- calcBinnedMotifEnrR(seqs = peak_seqs,bins =as.factor(input_bins) ,pwmL =pwms,verbose = TRUE)
monalisa_res
}

filter_motif_data <- function(monalisa_out,log10_p_thresh = 4,type = 'negLog10Padj'){
monalisa_filter <- apply(assay(monalisa_out, type), 1, 
             function(x) max(abs(x), 0, na.rm = TRUE)) >log10_p_thresh 
filtered_data <- monalisa_out[monalisa_filter,]

}

create_summary_data <- function(monalisa_output,type = 'negLog10Padj'){
motif_gene_mapping <- annotate_motifs(monalisa_output)
motif_summary_dat <- assay(monalisa_output,type) %>%
     data.frame() %>%
    rownames_to_column('motif') %>% 
    left_join(motif_gene_mapping,by = c('motif' = 'motif.id')) %>%
    distinct()
motif_summary_dat
}



pull_motifs_monalisa_output <- function(monalisa_output,motif_list){
selected_TF_indexes <- data.frame(motif = rowData(monalisa_output)$motif.name) %>%
    rownames_to_column('motif_id') %>%
    mutate(index = row_number()) %>%
    filter(toupper(motif) %in% motif_list | motif_id %in% motif_list) %>%
    pull(index)
selected_TFs <- monalisa_output[selected_TF_indexes,]
selected_TFs
}


# plots distributions for GC content and length for different bins 
# for an input dataframe
plot_bin_diagnostics <- function(input_df,bin_name,plot_name){
require(BSgenome.Hsapiens.UCSC.hg38)
pdf(plot_name)

input_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38,input_df %>% makeGRangesFromDataFrame())
input_bins <- input_df %>% pull({bin_name})  
plotBinDiagnostics(input_seqs,as.factor(input_bins),aspect = c("length"))
plotBinDiagnostics(input_seqs,as.factor(input_bins),aspect = c("GCfrac"))
dev.off()
}

# function that adjusts the length of granges by the median length
adjust_length <- function(input_granges){
adjusted <- trim(resize(input_granges, width = median(width(input_granges)), fix = "center"))
adjusted
}



annotate_motifs <- function(monalisa_output,mirror = "uswest"){
require(gprofiler2)
require(biomaRt)
#ensembl <-useMart("ensembl",mirror = mirror)
ensembl <- useEnsembl(biomart = 'ensembl',dataset = 'hsapiens_gene_ensembl',mirror = mirror)
#esemblist <- as.data.frame(listDatasets(ensembl))
#ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
t2g <- getBM(attributes=c('ensembl_gene_id','chromosome_name'), mart = ensembl)


motif_data <- rowData(monalisa_output) %>%
        data.frame() %>%
        dplyr::select(motif.id,motif.name)
gene_conversion_table <- gconvert(motif_data$motif.name) %>%
        dplyr::select(input,target,name) %>%
        left_join(motif_data,by = c('input' = 'motif.name')) %>%
        left_join(t2g,by = c( 'target'= 'ensembl_gene_id') )
gene_conversion_table
}

GetGC <- function(gr,bsgenome){
    seqs <- BSgenome::getSeq(bsgenome, gr)
    GC_vec <-  as.numeric(Biostrings::letterFrequency(x = seqs, letters = "GC", as.prob = TRUE))
    output <- gr %>%
        mutate(GC_content =GC_vec )
    output
}


plot_adjustment_comparison <- function(unadjusted_res,adjusted_res,type = 'negLog10Padj',plot = TRUE){
unadjusted_enrichment <- assay(unadjusted_res,type) %>%
    data.frame() %>%
    rownames_to_column('motifs') %>%
    pivot_longer(!motifs) %>%
    mutate(type = 'unadjusted')

pval_adjust_vs_unadjusted <- assay(adjusted_res,type) %>%
    data.frame() %>%
    rownames_to_column('motifs') %>%
    pivot_longer(!motifs) %>%
    mutate(type = 'adjusted') %>%
    bind_rows(unadjusted_enrichment) %>%
    pivot_wider(names_from = type,values_from = value)  
if (plot == TRUE) {
plot_out <-  pval_adjust_vs_unadjusted %>%
    #mutate(bin_type = case_when(name == 'FALSE.' ~ 'non-promoters',TRUE~'Promoters')) %>% 
    ggplot(aes(x = unadjusted,y = adjusted)) + 
        geom_point() +
        geom_hline(yintercept = 1.3,linetype = 'dashed',alpha = 0.5) + 
        geom_vline(xintercept = 1.3,linetype = 'dashed',alpha = 0.5) +
        xlab('P-values unadjusted bins') +
        ylab('P-values adjusted bins') + 
        theme_classic() + 
        facet_wrap(~name,nrow = 1) + 
        theme(strip.background = element_blank(),strip.text = element_text(size = 15))
plot_out
    }
else {
pval_adjust_vs_unadjusted
    }
}


