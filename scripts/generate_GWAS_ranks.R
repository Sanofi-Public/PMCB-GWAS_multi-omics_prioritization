#' Main function to generate the rankings from GWAS
#' #' It summarizes evidences for each gene and GWAS files considering
#' variant annotations (coding variants, ABC maps) and MR/Colocalization
#' It also provides locus to gene evidence rankings from "no evidence" to 
#' "very high" evidence, and flags genes with the highest evidence at a given
#' Locus. Each GWAS is treated separately.
#'

#' @param annotated_variants Variants annotated with VEP. Only significant variants should be included.
#' @param mr_coloc_results Results from mendelian randomization and colocalization analyses
#' @param abc_file Path to file containing ABC predictions
#' @param variant_pval_magnitude Only variants within this P-value magnitude of the lead variant will be considered
#' @param hla_start HLA region start (default=28477797)
#' @param hla_end HLA region start (default=33448354)
#'
#' @return A summary data table linking each Locus to genes
#' @export
#'
generate_GWAS_ranks = function(annotated_variants,
                            mr_coloc_results,
                            abc_file,
                            pathogenicity_file = NULL,
                            variant_pval_magnitude = 3, 
                            hla_start = 28477797, 
                            hla_end=33448354){
  library(plyr)
  library(data.table)
  library(GenomicRanges)
  library(dplyr)
  library(magrittr)
  
  results = vep_pathogenic_abc(annotated_variants, abc_file, pathogenicity_file)
  
  # Make sure all are complete
  results = results %>% mutate(
    is_ABC_TargetGene = ifelse(is_ABC_TargetGene=="" | is.na(is_ABC_TargetGene ), FALSE, is_ABC_TargetGene),
    is_coding = ifelse(is_coding=="" | is.na(is_coding ), FALSE, is_coding),
    is_pathogenic = ifelse(is_pathogenic=="" | is.na(is_pathogenic), FALSE, is_pathogenic)
  )
  
  # keep only protein coding genes and flag HLA genes
  summary = results %>% filter(BIOTYPE=="protein_coding") %>% mutate(HLA=ifelse(Chromosome==6 & Position > !!hla_start & Position< !!hla_end, TRUE, FALSE)) %>%
    dplyr::select(Gene, Symbol, GWAS_study_id, Phenotype, any_of("Sample_size"), any_of("N_cases"), HLA)
  summary = summary %>% distinct(Gene, GWAS_study_id, .keep_all=TRUE) 
  
  # add top SNP information to summary
  results.summary = results %>% dplyr::select(GWAS_study_id, Gene, Symbol, Top_SNP=SNP, Top_SNP_chr=Chromosome, Top_SNP_pos=Position, Top_SNP_beta=Beta, Top_SNP_se=SE, Top_SNP_eaf=EAF, Top_SNP_pval=Pvalue, Top_SNP_effect_allele=Alt, Top_SNP_other_allele=Ref, Distance)
  results.summary = results.summary %>% arrange(Top_SNP_pval) 
  results.summary = results.summary %>% mutate(Top_SNP_chr = as.character(Top_SNP_chr))
  results.summary = results.summary %>% distinct(Gene, GWAS_study_id, .keep_all=TRUE) 
  summary = summary %>% left_join(results.summary) %>% unique
  summary = summary %>% filter(!is.na(Gene))
  
  # Identify loci (aggregate variants that overlap within ~250kb)
  summary = annotate_loci_in_summary(summary)
  
  # flag HLA if any gene in the Locus overlaps the region
  hla = summary %>% group_by(GWAS_study_id, Locus) %>% dplyr::summarize(HLA = any(HLA, na.rm=TRUE))
  summary = summary %>% dplyr::select(-any_of("HLA")) %>% left_join(hla)
  
  # annotate the lead variant at each loci
  leads = summary %>% group_by(Locus, GWAS_study_id) %>% reframe(Top_SNP = Top_SNP[which(Top_SNP_pval == min(Top_SNP_pval, na.rm = TRUE))]) %>% mutate(Lead_SNP = TRUE) %>% unique
  summary = summary %>% left_join(leads) %>% replace_na(list(Lead_SNP=FALSE))
  
  # Identify the nearest gene to the lead SNP in each loci
  summary = annotate_nearest_gene(summary)
  
  # Prepare summary from annotations and MR/coloc results
  summary = summarize_features(summary, results, mr_coloc_results, variant_pval_magnitude = variant_pval_magnitude)
  
  # Rank each genes at each Locus
  summary = rank_summary(summary)
  
  # Find genes with highest scores
  summary = find_highest_rank_in_loci(summary)
  
  return(summary)
}

