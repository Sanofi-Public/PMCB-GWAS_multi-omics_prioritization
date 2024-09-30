
#' Rank genes at each Locus based on different rules
#' Highest evidence for Locus with lead coding variants or
#' ABC map support in addition to eQTL support
#'
#' @param summary summary data table with loci annotated
#'
#' @return summary data table annotated with evidence "ranks"
#' @export
#'
rank_summary = function(summary) {
  # Evidence - all
  summary = summary %>% mutate(
    Evidence_rank = case_when(
      (N_coloc_80 > 2 &
         Lead_SNP_ABCMax_linked) | Lead_Pathogenic ~ "VERY HIGH",
      (N_coloc_80 > 2 &
         Significant_SNP_ABCMax_linked) |
        (Coloc_H4_80 & Lead_SNP_ABCMax_linked)
      | Pathogenic | Lead_Lof_Missense ~ "HIGH",
      MR_pQTL |
        Coloc_H4_80
      | Lead_SNP_ABCMax_linked | Lead_Lof_Missense ~ "MODERATE" ,
      Significant_SNP_ABCMax_linked |
        Lead_SNP_ABC_linked | Nearest_gene_by_locus  ~ "WEAK",
      MR |
        Significant_SNP_ABC_linked ~ "VERY WEAK",
      .default = "NO EVIDENCE"
    )
  )
  
  summary = summary %>% mutate(
    Evidence_rank_numeric = recode(
      Evidence_rank,
      "VERY HIGH" = 5,
      "HIGH" = 4,
      "MODERATE" = 3,
      "WEAK" = 2,
      "VERY WEAK" = 1,
      "NO EVIDENCE" = 0
    )
  )
  
  summary %>% unique
  
}
