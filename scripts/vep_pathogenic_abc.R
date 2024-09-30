#' Extract annotations from VEP and external sources
#'
#' @param VEP_Results Data table with VEP results
#' @param abc_file ABC maps from Nasser, Nature 2021 . https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz
#' @param pathogenic_file Missense pathogenicity predictions (optional)
#'
#' @return Annotated results
#' @export
#'
vep_pathogenic_abc <- function(VEP_Results,
                               abc_file,
                               pathogenicity_file = NULL) {
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
  results = VEP_Results
  
  
  # read and process abc -----------------------------------------------------------------------
  abc = fread(abc_file)
  abc = abc %>% dplyr::select(chr, start, end, TargetGene, ABC.Score) %>%
    dplyr::mutate(chr = gsub("chr", "", chr)) %>%
    dplyr::mutate(ABC.Score = round(ABC.Score, 4))
  abc = makeGRangesFromDataFrame(abc, keep.extra.columns = TRUE)
  
  snps = results %>%
    dplyr::select(SNP,
                  Chromosome,
                  Position) %>%
    unique
  snps = makeGRangesFromDataFrame(
    snps,
    seqnames.field = "Chromosome",
    start.field = "Position",
    end.field = "Position",
    keep.extra.columns = TRUE
  )
  abc = mergeByOverlaps(snps, abc)
  abc = as.data.table(abc)
  abc_max = abc %>% group_by(SNP) %>% reframe(TargetGene = TargetGene[which(ABC.Score == max(ABC.Score, na.rm = TRUE))])
  abc = abc %>% dplyr::select(SNP, Symbol = TargetGene) %>% mutate(is_ABC_TargetGene = TRUE) %>% unique
  abc_max = abc_max %>% dplyr::select(SNP, Symbol = TargetGene) %>% mutate(is_ABCMax_TargetGene = TRUE) %>% unique
  
  # update results
  results = results %>%
    dplyr::rename(Symbol = SYMBOL) %>%
    left_join(abc) %>% left_join(abc_max)
  
  # Coding variants subsets
  # From VEP results
  results = results %>% mutate(is_coding =
                                 BIOTYPE == "protein_coding" &
                                 ((IMPACT == "MODERATE" &
                                     (
                                       grepl("missense", Consequence)
                                     )) | IMPACT == "HIGH"))
  
  if (!is.null(pathogenicity_file)) {
    pathogenicity = fread(pathogenicity_file)
    pathogenicity = pathogenicity %>%
      mutate(transcript_id = gsub("\\..*", "", transcript_id)) %>%
      mutate(Amino_acids = paste0(AA1, "/", AA2))
    pathogenicity2 = pathogenicity %>%  # in case its reversed
      mutate(Amino_acids = paste0(AA2, "/", AA1))
    pathogenicity = rbind(pathogenicity, pathogenicity2)
    pathogenicity = pathogenicity %>% distinct(Feature = transcript_id,
                                               Amino_acids,
                                               Protein_position = position) %>% mutate(is_pathogenic = TRUE)
    results = results %>% mutate(Protein_position = as.numeric(Protein_position)) %>%
      left_join(pathogenicity) %>% replace_na(list(is_pathogenic = FALSE))
  } else{
    results = results %>% mutate(is_pathogenic = FALSE)
  }
  
  ## Consider high impact variants
  results = results %>% mutate(is_pathogenic = IMPACT == "HIGH" |
                                 is_pathogenic)
  
  results = results %>% mutate(Distance = as.numeric(DISTANCE)) %>%
    mutate(Distance = ifelse(is.na(Distance) &
                               !is.na(Gene), 0, Distance))
  
  return(results)
}

