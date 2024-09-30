#' Generate a list of loci by aggregating variants in proximity 
#'
#' @param summary summary data table 
#' @param range genomic distance to consider (default=250000). Variants overlaping by this distance will be aggregated within the same Locus.
#'
#' @return annotated summary table
#' 
annotate_loci_in_summary = function(summary, range = 250000){
  
  loci = summary %>% group_by(GWAS_study_id) %>% find_top_genes(range=range)
  summary = summary %>% left_join(loci) 
  summary
  
}


#' Find top genes per locus
#'
#' @param study_results summary results for one GWAS
#' @param range genomic distance to consider (default=250000). Variants overlaping by this distance will be aggregated within the same Locus.
#'
#' @return top genes per locus
#'
#' @examples
find_top_genes = function(study_results, range = 250000){
  GWAS_study_id= unique(study_results$GWAS_study_id)
  sub_gr = makeGRangesFromDataFrame(study_results, seqnames.field = "Top_SNP_chr", start.field = "Top_SNP_pos", end.field = "Top_SNP_pos", keep.extra.columns = TRUE)+ range
  sub_reduced = GenomicRanges::reduce(sub_gr)
  top_genes = do.call(rbind, lapply(split(sub_reduced, as.factor(sub_reduced)), function(range){
    Locus = subsetByOverlaps(sub_gr, range)
    Chromosome = as.character(unique(seqnames(Locus)))
    start = min(start(Locus), na.rm=T)
    end = max(end(Locus), na.rm=T)
    data.table(GWAS_study_id= GWAS_study_id, Gene = Locus$Gene, Locus = paste(Chromosome, start, end, sep="_")) %>% unique
  }))
  top_genes
}