#'  Find highest ranking genes based on Evidence_rank
#'
#' @param summary a ranked summary table
#'
#' @return a summary table annotated with highest ranking genes at each Locus.
#' @export
#'
#' @examples
find_highest_rank_in_loci = function(summary) {
  ranked = summary %>% group_by(GWAS_study_id, Locus) %>% dplyr::reframe(Gene = unique(Gene[which(Evidence_rank_numeric ==
                                                                                                    max(Evidence_rank_numeric, na.rm = T))])) %>% ungroup
  ranked = ranked %>% mutate(Highest_ranked_gene_in_locus = TRUE)
  summary = summary %>% left_join(ranked) %>% mutate(Highest_ranked_gene_in_locus = ifelse(
    is.na(Highest_ranked_gene_in_locus),
    FALSE,
    Highest_ranked_gene_in_locus
  ))
  
  # ranked,but prioritized nearest
  ranked_near = summary %>% group_by(GWAS_study_id, Locus) %>% dplyr::reframe(Gene = if (# If there is more than one prioritized gene, and at least one gene is categorized as nearest gene, prioritize the latter,
    # unless there are high scoring genes in which case we keep all of them
    length(unique(Gene[which(Evidence_rank_numeric == max(Evidence_rank_numeric, na.rm =
                                                          TRUE))])) > 1 &
    length(unique(Gene[which(
      Evidence_rank_numeric == max(Evidence_rank_numeric, na.rm =
                                   TRUE) &
      Nearest_gene_by_locus
    )])) >= 1) {
    unique(Gene[which(
      Evidence_rank_numeric == max(Evidence_rank_numeric, na.rm =
                                     TRUE) &
        Nearest_gene_by_locus
    )])
  } else{
    # If there is only one prioritized gene or the lead variant is intergenic, we keep all prioritized genes
    unique(Gene[which(Evidence_rank_numeric == max(Evidence_rank_numeric, na.rm =
                                                     TRUE))])
  }) %>% ungroup
  ranked_near = ranked_near %>% mutate(Highest_ranked_gene_in_locus_plus_nearest =
                                         TRUE)
  summary = summary %>% left_join(ranked_near) %>% replace_na(list(Highest_ranked_gene_in_locus_plus_nearest = FALSE))
  
  return(summary)
  
}