#' Add features to summary based on annotations and MR/COLOC
#'
#' @param summary a summary data table
#' @param results annotated results data table
#' @param mr_coloc_res MR/COLOC results
#' @param variant_pval_magnitude threshold to use to consider annotations
#' 
#' @return a summary table with additional features
#' @export
#'
summarize_features = function(summary, results, mr_coloc_res, variant_pval_magnitude = 3){
  
  # add Locus info to results
  loci = summary %>% distinct(GWAS_study_id, Gene, Locus)
  
  results = results %>% inner_join(loci)
  summary_loci_leadSNPs = results %>% group_by(GWAS_study_id, Locus) %>% filter(Pvalue == min(Pvalue, na.rm=T)) %>% ungroup() %>% distinct(GWAS_study_id, Locus, Top_SNP=SNP) %>% mutate(Lead_SNP=TRUE) 
  summary = summary %>% left_join(summary_loci_leadSNPs) %>% replace_na(list(Lead_SNP=FALSE))
  results = results %>% left_join(summary_loci_leadSNPs, by=c("GWAS_study_id", "Locus", "SNP"="Top_SNP")) %>% replace_na(list(Lead_SNP=FALSE))
  
  # lead SNP extended
  lead_snps_extended = results %>% group_by(GWAS_study_id, Locus) %>% filter(
    -log10(Pvalue) >= -log10(min(Pvalue, na.rm=T))-!!variant_pval_magnitude
  ) %>% distinct(GWAS_study_id, Locus, Top_SNP=SNP) %>% mutate(Lead_SNP_extended=TRUE)
  results = results %>% left_join(lead_snps_extended, by=c("GWAS_study_id", "Locus", "SNP"="Top_SNP")) %>% replace_na(list(Lead_SNP_extended=FALSE))
  
  results = results  %>% 
    mutate(HIGH_IMPACT = IMPACT == "HIGH")
  
  # Annotate by Gene
  # TRUE if the condition is met at least once
  summarized_annotations = results %>% group_by(GWAS_study_id, Gene) %>%
    dplyr::summarize(
      Lof_Missense = any( is_coding & Lead_SNP_extended, na.rm = TRUE),
      Lead_Lof_Missense = any((is_coding) & Lead_SNP, na.rm=TRUE),
      Pathogenic = any(is_pathogenic & Lead_SNP_extended, na.rm=TRUE),
      Lead_Pathogenic = any(is_pathogenic & Lead_SNP, na.rm=TRUE),
      Significant_SNP_ABC_linked = any(is_ABC_TargetGene & Lead_SNP_extended, na.rm=TRUE),
      Lead_SNP_ABC_linked = any(is_ABC_TargetGene & Lead_SNP, na.rm=TRUE),
      Significant_SNP_ABCMax_linked = any(is_ABCMax_TargetGene & Lead_SNP_extended, na.rm=TRUE),
      Lead_SNP_ABCMax_linked = any(is_ABCMax_TargetGene & Lead_SNP, na.rm=TRUE)
    ) %>% ungroup()
  
  summary = summary %>% left_join(summarized_annotations) %>% distinct()
  
  # Add MR/Coloc information
  summarized_mr_coloc = mr_coloc_res %>% group_by(GWAS_study_id, Gene) %>% dplyr::summarise(
    MR = any(MR_qval  < 0.05, na.rm=T),
    MR_pQTL = any(MR_qval  < 0.05 & grepl("pQTL", Exposure, ignore.case=TRUE), na.rm=TRUE),
    MR_at_least_3_instruments = any(MR_qval  < 0.05 & MR_n_snp >=3, na.rm=TRUE),
    Coloc_H4_30 = any(MR_qval  < 0.05 & COLOC_PP_H4 >=0.3, na.rm=TRUE),
    Coloc_H4_80 = any(MR_qval  < 0.05 & COLOC_PP_H4 >=0.8, na.rm=TRUE),
    Coloc_H4_95 = any(MR_qval  < 0.05 & COLOC_PP_H4 >=0.95, na.rm=TRUE),
    Max_coloc_PPH4 = ifelse(!all(is.na(COLOC_PP_H4)), max(COLOC_PP_H4, na.rm = T), NA),
    N_coloc_80 = length(unique(Exposure_id[which(COLOC_PP_H4 > 0.80)])),
    N_coloc_95 = length(unique(Exposure_id[which(COLOC_PP_H4 > 0.95)]))
  )
  
  summary = summary %>% left_join(summarized_mr_coloc) %>% distinct()
  
  summary = summary %>% replace_na(
    list(
      MR = FALSE,
      Lof_Missense = FALSE,
      Lead_Lof_Missense = FALSE,
      Pathogenic = FALSE,
      Lead_Pathogenic = FALSE,
      Significant_SNP_ABC_linked = FALSE,
      Lead_SNP_ABC_linked = FALSE,
      Significant_SNP_ABCMax_linked = FALSE,
      Lead_SNP_ABCMax_linked = FALSE,
      MR_pQTL = FALSE,
      MR_at_least_3_instruments = FALSE,
      Coloc_H4_30 = FALSE,
      Coloc_H4_80 = FALSE,
      Coloc_H4_95 = FALSE,
      N_coloc_80 = 0,
      N_coloc_95 = 0
    )
  ) 
  
  return(summary)
  
}