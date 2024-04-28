# Takes in an ABC TSV file and filters on ABC score to high confidence 
# interactions and regions that are not promoter of gene. Then filters 
# to genes that exist in the dataframe of interest, so genes must be supplied 
# in input_granges, defaults to geneSymbol and then annotates which peaks are ABC regions.
# Currently drops the target genes of ABC regions, but that could be implemented back in to create 
# a column that contains the target genes as a list
annotate_ABC <- function(input_granges,ABC_path,ABC_thresh = 0.02,gene_column = 'geneSymbol'){
require(plyranges)
require(tidyverse)
gene_list <- input_granges %>% data.frame() %>% pull({gene_column})
ABC_data_filtered <- fread(ABC_path) %>% 
    filter(ABC.Score > ABC_thresh & isSelfPromoter == FALSE) %>%
    filter(TargetGene %in% gene_list) %>%
    dplyr::select(chr:end,CellType,-TargetGene) %>%
    distinct() %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE)

ABC_annotated <- input_granges %>%
    join_overlap_left(ABC_data_filtered) %>%
    data.frame() %>%
    mutate(ABC_peak = case_when(is.na(CellType) ~ FALSE,TRUE~TRUE)) %>%
    distinct() %>%
    pivot_wider(names_from = CellType,values_from = ABC_peak,values_fill = FALSE) %>% 
    dplyr::select(-`NA`) %>%
    distinct() %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE)
ABC_annotated
}




classify_peaks <- function(input_atac_path = '/oak/stanford/groups/smontgom/epadhi/ts21_motif_analysis/DownSyndrome_HSC_PeakGeneSets/pb_de_atac.hsc.txt',nominal_thresh = 0.05,adj_thresh = 0.05){
require(GenomicRanges)
require(data.table)
require(tidyverse)
all_peaks <- fread(input_atac_path) %>%
    separate(names,into = c('chr','start','end')) %>%
    mutate(direction = case_when(sign(logFC) == 1 ~ 'ts21/up', TRUE~'disomy/down')) %>%
    mutate(adj_sig = case_when(adj.P.Val < adj_thresh ~'sig',TRUE~'nonsig'),nominal_sig = case_when(P.Value < nominal_thresh ~'sig',TRUE~'nonsig')) %>% 
    dplyr::rename_with(.cols = c(1,2,3,4,5,6,10,11,12),~paste0(.,'_ATAC')) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) 
all_peaks
}

classify_expression <- function(input_rna_path = '/oak/stanford/groups/smontgom/epadhi/ts21_motif_analysis/DownSyndrome_HSC_PeakGeneSets/Liver.HSCs_MPPs.sample.txt'){
require(GenomicRanges)
require(data.table)
require(tidyverse)
expression_DE <- fread(input_rna_path) %>% 
    mutate(direction = case_when(sign(logFC) == 1 ~ 'ts21/up', TRUE~'disomy/down')) %>%
    mutate(adj_sig = case_when(adj.P.Val < .05 ~'sig',TRUE~'nonsig'),nominal_sig = case_when(P.Value < .05 ~'sig',TRUE~'nonsig')) %>% 
    dplyr::select(logFC,AveExpr,adj.P.Val,P.Value,names,direction,adj_sig,nominal_sig) %>%
    dplyr::rename_with(.cols = c(1:4,6:8),~paste0(.,'_expr')) 
expression_DE
}

load_clustered_motifs <- function(transfac_file = '/home/epadhi/meme/radial_trees/JASPAR_2022_matrix_clustering_vertebrates_CORE_cluster_root_motifs.tf',annotations = '/home/epadhi/meme/radial_trees/clusters_motif_names.tab'){
require(TFBSTools)
require(universalmotif)
clustered_motifs <- read_transfac(transfac_file)%>%
    convert_motifs(.,class ='TFBSTools-PWMatrix') %>%
    do.call(PWMatrixList, .)
cluster_annotations <- fread(annotations,header = FALSE)
output <- list(motifs = clustered_motifs,annotations = cluster_annotations)
output
}

annotate_EPD_promoters_expression <- function(expression_df,EPD_promoters = '/oak/stanford/groups/smontgom/epadhi/annotations/EPDnew_promoters/hg38/EPDnew_v6_hg38.bed'){
require(data.table)
require(tidyverse)
require(plyranges)
EPD_promoter_annotations_DE_genes <- fread(EPD_promoters) %>% 
    separate(V4,into = c('geneSymbol','promoter_number')) %>%
    dplyr::select(-promoter_number) %>%
    left_join(expression_df,by = c('geneSymbol' = 'names')) %>%
    filter(!is.na(P.Value_expr)) %>%
    dplyr::rename('chr' =1,'start' =2,'end' = 3) %>%
    distinct() %>%
    #mutate(start = start - 200,end = end + 200) %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE) 
EPD_promoter_annotations_DE_genes
}

annotate_ATAC_peaks_EPD_promoters <- function(ATAC_df,EPD_promoters){
# overlap ATAC seq fragments with EPD promoters to annotate 
# which regions are functioning as promoters and creates  granges object
EPD_ATAC_overlap_expression_annotated <- ATAC_df %>%
    join_overlap_left(EPD_promoters) %>%
    mutate(promoter = case_when(is.na(geneSymbol) ~ FALSE,TRUE~TRUE)) %>% 
    dplyr::select(-V5,-V6,-V7,-V8) %>%
    data.frame() %>%
    distinct() %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE)
EPD_ATAC_overlap_expression_annotated

}


# runs go enrichment analysis and binds results on a column
run_GO_enrich <- function(input_df,gene_column = 'geneSymbol'){
require(enrichR)
dbs <- c( "GO_Cellular_Component_2023", "GO_Biological_Process_2023","GO_Molecular_Function_2023")
enrich_res <- enrichr(input_df %>% data.frame() %>% pull({gene_column}),dbs) %>%
    bind_rows(.id = 'column_label')
enrich_res
}

