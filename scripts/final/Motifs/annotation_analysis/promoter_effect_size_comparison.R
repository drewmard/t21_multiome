library(tidyverse)
library(data.table)
library(introdataviz)
library(SummarizedExperiment)
library(JASPAR2022)
library(TFBSTools)
library(motifmatchr)
source('../utils/data_processing.R')
sessionInfo()

JASPAR_all_motifs <- getMatrixSet(JASPAR2022)
pwms <- getMatrixSet(JASPAR2022,opts = list(matrixtype = "PWM",tax_group = "vertebrates"))

ABC_all_path <- '/oak/stanford/groups/smontgom/epadhi/ts21_motif_analysis/ABC_data/HSCsts21_EnhancerPredictionsAllPutative.tsv.gz'
all_peaks <- classify_peaks(adj_thresh = .1) 
expression_DE <- classify_expression() 
EPD_promoter_annotations_DE_genes <- annotate_EPD_promoters_expression(expression_DE) 
EPD_ATAC_overlap_expression_annotated <- annotate_ATAC_peaks_EPD_promoters(all_peaks,EPD_promoter_annotations_DE_genes) %>% 
    annotate_ABC(ABC_all_path)
   
ABC_data <- fread(ABC_all_path) %>%
    filter(ABC.Score > 0.02) %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE)

motif_metadata <- rowData(readRDS('../ts21_differential_motif_analysis/output_data/enhancer_ts21_vs_disomy_motif_length_adjust.RDS'))$motif.name %>%
    data.frame() %>%
    rownames_to_column('JASPAR_id') %>% 
    dplyr::rename('motif_name' = 2)
########ANALZYE NUMBER OF ABC ENHANCESR PER GENE 
plot_ABC_per_gene <- function(ABC_ATAC_peaks_joined){
non_ABC_genes <- expression_DE %>%
    filter(!names %in% ABC_list) %>%
    select(names) %>%
    distinct() 
ABC_filtered_data <- ABC_ATAC_peaks_joined %>%
    filter(isSelfPromoter == FALSE) %>%
    group_by(TargetGene) %>% 
    summarize(num_ABC = dplyr::n()) %>%
    #dplyr::count(TargetGene) %>%
    #dplyr::count(n) %>%
    rbind(c(0,length(non_ABC_genes$names)))  %>%
    ggplot(aes(x =n, y =nn))  +
        geom_col() + 
        theme_classic() + 
        xlab('Number of ABC enhancers per Gene') +
        ylab('Number of Genes') + 
        scale_y_continuous(expand = c(0,0),breaks = c(0,500,1000,1500,2000,2500,3000))
ABC_filtered_data
}

ABC_enhancers_per_gene <- ABC_data %>%
    filter(isSelfPromoter == FALSE) %>%
    join_overlap_left(all_peaks) %>%
    data.frame() %>% 
    inner_join(expression_DE,by = c('TargetGene' = 'names')) %>%
    filter(!is.na(logFC_expr) & !is.na(adj_sig_ATAC)) %>%
    dplyr::count(TargetGene) %>%
    dplyr::count(n) %>%
    rbind(c(0,length(non_ABC_genes$names)))  %>%
    ggplot(aes(x = n,y = nn)) + 
        geom_col() 


DA_ABC_enhancers_per_gene <-  ABC_data %>%
    filter(isSelfPromoter == FALSE) %>%
    join_overlap_left(all_peaks) %>%
    data.frame() %>% 
    inner_join(expression_DE,by = c('TargetGene' = 'names')) %>%
    filter(!is.na(logFC_expr) & !is.na(adj_sig_ATAC)) %>%
    group_by(TargetGene) %>%
    summarize(sig = sum(adj_sig_ATAC == 'sig')) %>%
    dplyr::count(sig) %>% 
    ggplot(aes(x = as.factor(sig),y = log10(n))) + 
        geom_col() +
        theme_classic() + 
        xlab('Number of DA ABC enhancers per Gene') +
        ylab('log10(Number of Genes)') + 
        scale_y_continuous(expand = c(0,0)) 
ABC_expression_joint <- ABC_data %>%
    join_overlap_left(all_peaks) %>%
    data.frame() %>%
    inner_join(expression_DE,by = c('TargetGene' = 'names')) %>%
    filter(!is.na(logFC_expr) & !is.na(adj_sig_ATAC)) 

# calculate the association between differential acessibility and how 
# likely those enhancers are to be connected to a DE gene 
enhancer_expression_enrichment <- ABC_expression_joint %>%
    filter(isSelfPromoter == FALSE) %>%
    xtabs(~ adj_sig_ATAC + adj_sig_expr,data = .) %>% 
    fisher.test() %>%
    broom::tidy() %>%
    mutate(type = 'Observing gene DE \n given ABC enhancer DA')

# calculate the association between observing a differentially expressed 
# gene and observing atleast one DA ABC enhancer
gene_expression_enhancer_enrichment <- ABC_expression_joint %>%
    filter(isSelfPromoter == FALSE) %>% 
    mutate(boolean_adj_sig_ATAC = case_when(adj_sig_ATAC == 'sig' ~ TRUE,TRUE~ FALSE),isDE  = case_when(adj_sig_expr == 'sig' ~ TRUE,TRUE~ FALSE)) %>%
    group_by(TargetGene) %>%
    summarize(hasDAABC = sum(boolean_adj_sig_ATAC == TRUE),isDE = isDE,logFC = logFC_expr) %>% 
    ungroup() %>%
    distinct()  %>%
    mutate(hasDAABC = case_when(hasDAABC > 0 ~ TRUE,TRUE ~ FALSE)) %>%
    xtabs(~hasDAABC + isDE ,data = .) %>% 
    fisher.test() %>%
    broom::tidy() %>%
    mutate(type = 'Observing >= 1 DA ABC \n enhancers given gene DE')

ABC_expression_joint %>%
    filter(isSelfPromoter == FALSE) %>% 
    mutate(boolean_adj_sig_ATAC = case_when(adj_sig_ATAC == 'sig' ~ TRUE,TRUE~ FALSE),isDE  = case_when(adj_sig_expr == 'sig' ~ TRUE,TRUE~ FALSE)) %>%
    group_by(TargetGene) %>%
    summarize(hasDAABC = sum(boolean_adj_sig_ATAC == TRUE),isDE = isDE,logFC = logFC_expr) %>% 
    ungroup() %>%
    distinct()  %>%
    mutate(hasDAABC = case_when(hasDAABC > 0 ~ TRUE,TRUE ~ FALSE)) %>%
    xtabs(~hasDAABC + isDE ,data = .) 

OR_plots <- gene_expression_enhancer_enrichment %>%
    bind_rows(enhancer_expression_enrichment) %>%
    ggplot(aes(y =estimate, x = type )) +
        geom_errorbar(aes(ymin = conf.low,ymax = conf.high )) +
        geom_point(size = 8,color = 'orange') +
        theme_classic() +
        ylab('Odds Ratio') + 
        geom_hline(yintercept = 1,linetype = 'dashed',color = 'red',alpha = 0.7)

ABC_plot <- ABC_enhancers_per_gene | DA_ABC_enhancers_per_gene | OR_plots
ggsave('ABC_enhancers_per_gene.pdf',plot = ABC_plot,width = 12)



###### COMPARE EFFECT SIZES OF PROMOTERS VS ABC
# proportion of EPD promoters that are differentially accessibile 
EPD_ATAC_overlap_expression_annotated %>%
    filter(promoter == TRUE ) %>%
    data.frame() %>%
    dplyr::count(adj_sig_ATAC) 

# Proportion of ABC enhancers that are differentially 
# accessible
EPD_ATAC_overlap_expression_annotated %>%
    filter(HSCsts21 == TRUE & promoter == FALSE ) %>%
    data.frame() %>%
    dplyr::count(adj_sig_ATAC) 


logFC_ATAC_ABC_vs_EPD <- EPD_ATAC_overlap_expression_annotated %>%
    filter(promoter == TRUE | HSCsts21 == TRUE) %>%
    filter(adj_sig_ATAC == 'sig') %>%
    mutate(promoter = case_when(promoter == TRUE ~ 'EPD Promoters',TRUE~'ABC Enhancers')) %>%
    data.frame() %>%
    mutate(direction_ATAC = case_when(str_detect(direction_ATAC,'ts21') ~ 'Ts21',TRUE~'Disomy')) %>% 
    ggplot(aes(x = direction_ATAC,y = abs(logFC_ATAC),fill = promoter)) +
        geom_split_violin(alpha = 0.85,trim = FALSE) +
        #geom_violin(show.legend = FALSE) + 
        geom_boxplot(width = .08,alpha = 0.6,show.legend =  FALSE) +
        scale_fill_carto_d(palette = 'Safe',labels = c('ABC Enhancer','EPD Promoter'),name = 'Type') +  
        ylab('|ATAC logFC|') + 
        xlab('ATAC Direction') + 
        theme_classic() + 
        theme(legend.position = c(0.15,0.93),axis.text = element_text(size = 18),axis.title = element_text(size = 18),legend.text=element_text(size=14))  
ggsave('ABC_vs_EPD_effect_sizes.pdf')


EPD_ATAC_overlap_expression_annotated %>%
    filter(promoter == TRUE | HSCsts21 == TRUE) %>%
    filter(adj_sig_ATAC == 'sig') %>% 
    data.frame() %>%
    do(broom::tidy(lm(abs(logFC_ATAC) ~ promoter + direction_ATAC + promoter:direction_ATAC,data = .)))


sig_count <- EPD_ATAC_overlap_expression_annotated %>%
    filter(promoter == TRUE | HSCsts21 == TRUE) %>%
    filter(adj_sig_ATAC == 'sig') %>%
    data.frame() %>%
    dplyr::count(promoter,direction_ATAC) %>%
    ggplot(aes(x = direction_ATAC,fill = promoter, y = n)) +
        geom_col(position = 'dodge',color = 'black',alpha = 0.85) + 
        scale_fill_carto_d(palette = 'Safe',labels = c('ABC Enhancer','EPD Promoter'),name ='Annotation') + 
        ylab('Number significant') + 
        scale_y_continuous(expand = c(0,0)) +
        xlab('Direction ATAC logFC') + 
        theme_classic()
output <- logFC_ATAC_ABC_vs_EPD + sig_count
ggsave('ATAC_dist.pdf',width = 10)

enhancer_ts21_vs_disomy_motif_enr <- readRDS('../ts21_differential_motif_analysis/output_data/enhancer_ts21_vs_disomy_motif_length_adjust.RDS')
GATA_motif_list <- rowData(enhancer_ts21_vs_disomy_motif_enr)[,c('motif.name','motif.id')] %>% 
    data.frame() %>% 
    filter(str_detect(toupper(motif.name),'GATA')) %>% 
    pull(motif.id)
GATA_PWMs <- getMatrixByID(JASPAR2022,GATA_motif_list) 


matches_all_pwms <- EPD_ATAC_overlap_expression_annotated %>%
    filter( HSCsts21 == TRUE | promoter == TRUE) %>%
    filter(adj_sig_ATAC == 'sig') %>%
    matchMotifs(pwms,.,genome = 'hg38',out = 'scores')

motf_hit_matrix <- assays(matches_all_pwms)$motifScores %>% 
    data.matrix() %>%
    data.frame() %>% 
    mutate(across(everything(),~case_when(. > 1 ~ TRUE,TRUE~FALSE)))

all_motif_hits <- EPD_ATAC_overlap_expression_annotated %>%
    filter( HSCsts21 == TRUE | promoter == TRUE) %>%
    filter(adj_sig_ATAC == 'sig') %>%
    data.frame() %>%
    bind_cols(motf_hit_matrix)

model_res <- all_motif_hits %>%
    pivot_longer(contains('MA')) %>%
    group_by(name) %>%
    do(broom::glance(lm(logFC_ATAC ~ value + promoter  + value:promoter,data =.))) 


ATAC_permutations <- all_motif_hits %>% data.frame() %>% permute(.,1,contains('MA'))
get_TF_results <- function(input_df){
model_res <- input_df %>% 
    pivot_longer(contains('MA')) %>%
    group_by(name) %>%
    do(broom::glance(lm(logFC_ATAC ~ value + promoter  + value:promoter,data =.)))
out <- data.frame(max_r2 = max(model_res$adj.r.squared))
out
}

library(rcartocolor)
permutaitons <- ATAC_permutations %>%pmap_dfr(~as.data.frame(.) %>% get_TF_results)
permutaitons %>%
    arrange(desc(max_r2))
TF_R2 <- model_res %>%
    mutate(padj = p.adjust(p.value,method = 'BH')) %>%
    filter(adj.r.squared >= .035) %>%
    left_join(motif_metadata,by = c('name' = 'JASPAR_id')) %>%
    filter(motif_name != 'CTCF') %>%
    ggplot(aes(y = reorder(motif_name,r.squared),x = r.squared,fill = r.squared)) + 
        geom_col(show.legend = FALSE) + 
        ylab('Motif') +
        scale_fill_carto_c(palette = 'Teal') + 
        xlab(bquote(''~R^2)) + 
        theme_classic() +  
        scale_x_continuous(expand = expansion(mult = c(0,.2))) + 
        theme(axis.text = element_text(size = 23),axis.title = element_text(size = 25)) + 
        coord_cartesian(clip = "off")
ggsave('R2_TF_model.pdf',height = 15,width = 8)

matches_GATA <- EPD_ATAC_overlap_expression_annotated %>%
    filter( HSCsts21 == TRUE | promoter == TRUE) %>%
    filter(adj_sig_ATAC == 'sig') %>%
    matchMotifs(GATA_PWMs,.,genome = 'hg38',out = 'scores')
gata_motif_hits <- assays(matches_GATA)$motifScores %>% 
    data.matrix() %>%
    data.frame() %>% 
    mutate(across(everything(),~case_when(. > 1 ~ TRUE,TRUE~FALSE)))  %>%
    mutate(GATA_containing = case_when(if_any(everything(),~. == TRUE) ~ TRUE,TRUE~FALSE)) %>%
    pull(GATA_containing)

GATA_annotated <- EPD_ATAC_overlap_expression_annotated %>%
    filter(HSCsts21 == TRUE | promoter == TRUE) %>%
    filter(adj_sig_ATAC == 'sig') %>%
    data.frame() %>%
    mutate(GATA_containing = gata_motif_hits) 

GATA_model <- lm(scale(logFC_ATAC) ~ GATA_containing + promoter + GATA_containing:promoter,data =GATA_annotated %>% data.frame())
GATA_model <- lm(abs(logFC_ATAC) ~ GATA_containing + promoter + GATA_containing:promoter,data =GATA_annotated %>% data.frame())
broom::tidy(lm(abs(logFC_ATAC) ~ GATA_containing*promoter*direction_ATAC,data =GATA_annotated %>% data.frame()))
broom::tidy(GATA_model)


custom_pal <- carto_pal(8,'RedOr')[c(2,6)] 
GATA_plot <- GATA_annotated %>%
    data.frame() %>%
    mutate(annotation = case_when(promoter == TRUE ~ 'EPD Promoters',TRUE ~ 'ABC Enhancers')) %>%
    mutate(direction_ATAC = case_when(str_detect(direction_ATAC,'ts21') ~ 'Ts21',TRUE~'Disomy')) %>% 
    ggplot(aes(x = direction_ATAC,y = abs(logFC_ATAC),fill = GATA_containing)) +
        geom_split_violin(alpha = 0.75,trim = FALSE) +
        geom_boxplot(width = .1,alpha = 0.6,show.legend =  FALSE) +
        scale_fill_manual(values=custom_pal,name = 'GATA Containing') +
        ylab('|ATAC logFC|') + 
        xlab('Direction ATAC logFC') + 
        theme_classic() +
        facet_wrap(~annotation,nrow = 1) + 
        theme(legend.position = c(0.12,0.92),strip.background = element_blank(),legend.title = element_text(size=12),strip.text = element_text(size = 16),axis.text = element_text(size = 16),
              axis.title = element_text(size = 18),legend.text = element_text(size = 12)) 
ggsave('GATA_effect_size_comparison.pdf',width = 8)
#GATA_model <- lm(scale(logFC_ATAC) ~ GATA_containing + promoter + GATA_containing:promoter,data =GATA_annotated %>% data.frame())
summary(lm(abs(logFC_ATAC) ~  GATA_containing + promoter + GATA_containing:promoter ,data =GATA_annotated %>% data.frame()))



filtered_ABC <- ABC_data %>%
    filter(ABC.Score > .02) %>%
    makeGRangesFromDataFrame(keep.extra.columns =  TRUE)
GATA_ABC_target_genes <- filtered_ABC %>%
    join_overlap_left(GATA_annotated) %>%
    data.frame() %>%
    filter(GATA_containing == TRUE & !is.na(logFC_ATAC) & adj_sig_ATAC == 'sig' & promoter == FALSE ) 
GATA_up_ATAC_genes <- GATA_ABC_target_genes %>%
    filter(sign(logFC_ATAC) == 1 ) %>%
    distinct(TargetGene) %>%
    select(TargetGene) %>% 
    run_GO_enrich(gene_column ='TargetGene') %>%
    mutate(type = 'Trisomy 21')

GATA_down_ATAC_genes <- GATA_ABC_target_genes %>%
    filter(sign(logFC_ATAC) == -1 ) %>%
    distinct(TargetGene) %>%
    select(TargetGene) %>% 
    run_GO_enrich(gene_column ='TargetGene') %>%
    mutate(type = 'Disomy')

   
 
go_terms_up <- GATA_up_ATAC_genes %>%
    bind_rows(GATA_down_ATAC_genes) %>%
    filter(Adjusted.P.value < .1  ) %>%
    separate(Term,into = c('Pathway','GO'),sep = '\\(') %>%
    filter(Odds.Ratio > 8 & type == 'Trisomy 21' | Odds.Ratio > 23 & type == 'Disomy') %>%
    ggplot(aes(y = reorder(Pathway,Odds.Ratio),x = Odds.Ratio)) +
        geom_col(fill = custom_pal[2],color = 'black',alpha = 0.8) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_discrete(labels = scales::label_wrap(27)) +
        theme_classic() + 
        ylab('GO Pathway') + 
        xlab('Pathway Enrichment (OR) \n GATA regulated genes') + 
        facet_wrap(~type,nrow = 2,scales = 'free') + 
        theme(strip.background = element_blank(),strip.text = element_text(size = 15))
new_out <- (output / GATA_plot)
#out_again <- new_out  | go_terms_up 
ggsave('GATA_dist.pdf',height = 10,width = 10)


