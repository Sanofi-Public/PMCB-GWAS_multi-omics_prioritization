#' Function to annotate nearest genes to lead variants at each Locus
#'
#' @param summary summary data table
#'
#' @return summary data table with "nearest" genes annotated
#' @export
#'
annotate_nearest_gene = function(summary_df) {
  summary_df = summary_df %>% mutate(Top_SNP_chr = as.character(Top_SNP_chr))
  
  nearest = summary_df %>% group_by(Locus, GWAS_study_id) %>%
    filter(Distance == min(Distance, na.rm = T)) %>% mutate(Lead_SNP = TRUE) %>%
    ungroup() %>% distinct(Locus, GWAS_study_id, Gene, Top_SNP) %>%
    mutate(Nearest_gene_by_locus = TRUE)
  summary_df = summary_df %>% left_join(nearest)
  summary_df = summary_df %>% replace_na(list(Nearest_gene_by_locus = FALSE))
  summary_df
}