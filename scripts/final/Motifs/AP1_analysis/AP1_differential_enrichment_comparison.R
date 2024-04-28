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

create_base_motif_plot <- function(input_ggplot){
output <- input_ggplot + 
            geom_point(size = 3,alpha = 0.6) +
            geom_hline(yintercept = 1.3,linetype = 'dashed',alpha = 0.6) + 
            geom_vline(xintercept = 1.3,linetype = 'dashed',alpha = 0.6) +
            theme_classic() +
            xlab('-log10(Promoter p-val)') + ylab('-log10(Enhancer p-val)') +
            rcartocolor::scale_colour_carto_d(palette = 'Temps') + 
            theme(legend.position = c(0.85,0.92)) 
output 
}




######### LOAD ENRICHMENT DATA
#enhancer_ts21_vs_disomy_motif_enr <- readRDS('../ts21_differential_motif_analysis/output_data/enhancer_ts21_vs_disomy_motif_length_adjust.RDS')
ABC_enhancer_ts21_vs_disomy_motif_enr <- readRDS('../ts21_differential_motif_analysis/output_data/ABC_ts21_vs_disomy_motif_length_adjust.RDS')
ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig <- readRDS('../ts21_differential_motif_analysis/output_data/ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig.RDS')
ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig <- readRDS('../ts21_differential_motif_analysis/output_data/ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig.RDS')
ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background <- readRDS('../ts21_differential_motif_analysis/output_data/ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background.RDS')

AP1_motif_list <- rowData(ABC_enhancer_ts21_vs_disomy_motif_enr)$motif.name %>%
    data.frame() %>%
    filter(str_detect(.,'JUN|FOS|ATF')) %>%
    pull(.)

AP1_JASPAR_IDs <- rownames(pull_motifs_monalisa_output(ABC_enhancer_ts21_vs_disomy_motif_enr,AP1_motif_list)) %>%
    data.frame() %>%
    filter(. != 'MA1143.1') %>%
    pull(.)

all_AP1 <-getMatrixByID(JASPAR2022,ID=AP1_JASPAR_IDs)
all_ap1_cor <- get_motif_correlations(names(all_AP1),plot = FALSE)
all_ap1_clust <- hclust(dist(all_ap1_cor))

AP1_clusters <- cutree(all_ap1_clust,k = 2) %>%
        data.frame() %>% 
        rownames_to_column('motifs') %>%
        dplyr::rename('cluster' = 2) %>%
        mutate(AP1_motif = case_when(cluster == 2 ~ 'CRE Motif',TRUE ~ 'TRE Motif' ))

########
ABC_enhancer_AP1_data  <- pull_motifs_monalisa_output(ABC_enhancer_ts21_vs_disomy_motif_enr,AP1_motif_list)
#enhancer_AP1_data <- pull_motifs_monalisa_output(enhancer_ts21_vs_disomy_motif_enr,AP1_motif_list)
promoter_AP1_data_adj <- pull_motifs_monalisa_output(ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig,AP1_motif_list)
promoter_AP1_data_adj_background <- pull_motifs_monalisa_output(ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background,AP1_motif_list)
promoter_AP1_data_nom <- pull_motifs_monalisa_output(ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig,AP1_motif_list)


enrichments_promoter_adj <- plot_adjustment_comparison(promoter_AP1_data_adj,ABC_enhancer_AP1_data,plot = FALSE) %>%
                left_join(AP1_clusters,by ='motifs' ) %>%
                na.omit() %>% 
                ggplot(aes(x = unadjusted,y = adjusted,color = name,shape = AP1_motif)) %>% 
                    create_base_motif_plot() + 
                    labs(colour = "Enrichment Type",shape = 'AP-1 Motif Type')  
ggsave('plots/motif_comparisons/AP1_motif_comparison_adj_promoters.pdf',width = 5,height = 7)

enrichments_promoter_adj_background <- plot_adjustment_comparison(promoter_AP1_data_adj_background,ABC_enhancer_AP1_data,plot = FALSE) %>%
                left_join(AP1_clusters,by ='motifs' ) %>%
                separate(name,into = c('name','trash')) %>%
                mutate(name = case_when(name == 'ts21' ~ 'Ts21', TRUE ~ 'Disomy')) %>%  
                na.omit() %>% 
                ggplot(aes(x = unadjusted,y = adjusted,color = name,shape = AP1_motif)) %>% 
                    create_base_motif_plot() + 
                    labs(colour = "Enrichment Type",shape = 'AP-1 Motif Type') +
                    theme(legend.position = c(0.85,0.85)) 
ggsave('plots/motif_comparisons/AP1_motif_comparison_adj_promoters_background_seqs.pdf',width = 5,height = 6)


enrichments_promoter_nominal <- plot_adjustment_comparison(promoter_AP1_data_nom,enhancer_AP1_data,plot = FALSE) %>%
                left_join(AP1_clusters,by ='motifs' ) %>%
                na.omit() %>% 
                ggplot(aes(x = unadjusted,y = adjusted,color = name,shape = AP1_motif)) %>% 
                    create_base_motif_plot() + 
                    labs(colour = "Enrichment Type",shape = 'AP-1 Motif Type')  
ggsave('plots/motif_comparisons/AP1_motif_comparison_nominal_promoters.pdf',width = 5,height = 7)





motif_id_name_mapping <- create_summary_data(ts21_vs_disomy_EPDpromoters_ATAC_direction_nominal_sig) %>%
    dplyr::select(motif,name) %>%
    dplyr::rename('tf_name' =2)

enrichments_all_motifs <- plot_adjustment_comparison(ts21_vs_disomy_EPDpromoters_ATAC_direction_adj_sig_background,ABC_enhancer_ts21_vs_disomy_motif_enr,plot = FALSE) %>%
                filter(name != 'nonsig') %>%
                left_join(motif_id_name_mapping,by = c('motifs' = 'motif')) %>%
                separate(name,into = c('name','trash')) %>%
                mutate(name = case_when(name == 'ts21' ~ 'Ts21', TRUE ~ 'Disomy')) %>% 
                ggplot(aes(x = unadjusted,y = adjusted,color = name)) %>% 
                    create_base_motif_plot() + 
                    geom_text_repel(data =. %>% filter(unadjusted > 1.3 & adjusted > 100),aes(label = tf_name),show.legend = FALSE,box.padding = 0.5) + 
                    geom_text_repel(data =. %>% filter( adjusted > 10),aes(label = tf_name),show.legend = FALSE,box.padding = 0.5) + 
                    geom_text_repel(data =. %>% filter( unadjusted > 9),aes(label = tf_name),show.legend = FALSE,box.padding = 0.5,max.overlaps = 12) + 
                    xlab('-log10(Promoter p-val)') + ylab('-log10(Enhancer p-val)') +
                    labs(colour = "Enrichment Type")  
ggsave('all_motifs_comparison.pdf',width = 5,height = 6)
zoomed <- enrichments_all_motifs + ylim(0,100)
ggsave('all_motifs_comparison_zoomed.pdf',width = 5,height = 7,plot= zoomed)



