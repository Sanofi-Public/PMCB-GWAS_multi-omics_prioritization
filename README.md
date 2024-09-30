# Leveraging large-scale multi-omics to identify therapeutic targets from genome-wide association studies 

The code contained in this repository can be used to prioritize genes in GWAS loci using variant annotations, epigenomic data, and molecular quantitative trait loci (molQTL) as in Lessard *MedRxiv* 2023, doi [10.1101/2023.11.01.23297926](https://doi.org/10.1101/2023.11.01.23297926).

## Description

The script takes as input variants annotated using the [Variant Effect Predictor](https://ensembl.org/info/docs/tools/vep/index.html), results from mendelian randomization analyses using molQTL as exposure as well as colocalization analyses, activity-by-contact (ABC) maps, and optionally variant pathogenicity predictions. Each gene is annotated with these features and ranked using a rule-based approach (Lessard *MedRxiv* 2023, doi [10.1101/2023.11.01.23297926](https://doi.org/10.1101/2023.11.01.23297926)), prioritizing genes supported by both molQTL and maximum ABC scores, or the presence of associated coding variants.

## Input

The following inputs are required to run the script.

### 1. Input variants

The first required input is a dataframe that contains information about GWAS variants. The dataframe should only include variants below a certain p-value threshold ($p<5*10^{-6})$ or lower is recommended). Each row should contain the following information:

- GWAS_study_id: Unique GWAS ID. The input can contain multiple GWAS with different IDs.
- Phenotype: Phenotype tested in the GWAS
- SNP: Variant rsID
- Chromosome: Chromosome where the variant is located
- Position: Position of the variant
- Ref: Reference allele
- Alt: Alternative (effect) allele
- Beta: Variant effect
- SE: Standard error
- Pvalue: Association P-value
- EAF: Effect allele frequency
- **VEP annotations** (see below)

The dataframe should contain information from [Variant Effect Predictor (VEP)](https://ensembl.org/info/docs/tools/vep/index.html) (McLaren *Genome Biology*, 2016, doi [10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)), which has been run on the same GWAS with the options *--everything* and *--distance 250000* (or higher), for example:

```
perl vep -i in.vcf -o out.vcf --vcf --distance 250000 --cache --dir_cache /path/to/cache/directory --everything --offline --check_existing
```

The [parseCSQToGRanges](https://www.rdocumentation.org/packages/ensemblVEP/versions/1.12.0/topics/parseCSQToGRanges) function can be used to extract information from the VEP output:

```
# in R
out <- ensemblVEP::parseCSQToGRanges(VariantAnnotation::readVcf(outvcf, "hg19"), 
  VCFRowID = rownames(invcf)) %>% dplyr::as_tibble()
```

### 2. Mendelian randomization and colocalization results

The script requires results from mendelian randomization (MR) analyses using molQTL as exposures, in addition to colocalization results. Tools such as [gwasglue](https://github.com/MRCIEU/gwasglue) can be useful to generate the required inputs, and molQTL data can be retrieved from repositories like the [eQTL catalogue](https://www.ebi.ac.uk/eqtl/) (Kerimov *Nature Genetics*, 2021, doi [10.1038/s41588-021-00924-w](https://doi.org/10.1038/s41588-021-00924-w)).

The dataframe should minimally contain the following columns:

- GWAS_study_id: Unique GWAS ID (same as above)
- Exposure_id: Unique exposure ID. This corresponds to a unique molQTL dataset.  
- Gene: Ensembl gene ID corresponding to the tested molQTL.
- Exposure: Type of exposure (e.g. cis-eQTL, pQTL)
- MR_qval: MR adjusted p-value
- MR_n_snp: Number of variants included in the MR analysis
- COLOC_n_snp: Number of SNPs included in the colocalization analysis
- COLOC_PP_H4: Posterior probability of colocalization

### 3. Activity-by-contact (ABC) maps

ABC maps from e.g. [Nasser et al. 2021](https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz). The same format is expected.

### 4. Variant pathogenicity predictions (optional)

A file with genome-wide missense variant pathogenicity predictions can be used. These can be retrieved from different sources including [ProtVar](https://www.ebi.ac.uk/ProtVar/). **All variants in the file will be considered pathogenic**. The file should have the following columns:

- transcript_id: Ensembl transcript identifier (ENST*)
- AA1: Reference amino acid (One letter code)
- AA2: Alternative amino acid (One letter code)
- position: Protein amino acid position


## Running the code

To run the prioritization script:

```
# in R
library(data.table)

# source all scripts
scripts = list.files("./scripts/", pattern="*.R$",  full.names = TRUE)
sapply(scripts, source, .GlobalEnv)

# Read input data
annotated_variants = fread("/path/to/annotated/variants")
mr_coloc_results = fread("/path/to/MR/coloc/results")
abc_file_path = "/path/to/abc/file"
pathogenicity_file_path = "/path/to/pathogenicity/predictions" 
 
ranked_gwas = generate_GWAS_ranks(
  annotated_variants = annotated_variants, # VEP annotated variants 
  mr_coloc_results = mr_coloc_results, # Results from MR/coloc analyses 
  abc_file = abc_file_path, # Path to ABC file
  pathogenicity_file = pathogenicity_file, # Path to pathogenicity prediction file. Set to NULL is not used.
  variant_pval_magnitude = 5, # Only variants within this P-value magnitude from the lead variant will be considered.
  hla_start = 28477797, # HLA region start position (change if GWAS are on build hg38), 
  hla_end= 33448354 # HLA region end position (change if GWAS are on buld hg38)
)

```

## Output

The `generate_GWAS_ranks()` function outputs the final gene predictions. Each row contains information about a specific gene, for a specific GWAS. Several columns corresponding to the input features are reported. The *Evidence_rank* column reports the prioritization rank based on each features, from no evidence to very high evidence. The *Highest_ranked_gene_in_locus* column reports whether the gene is the gene with highest rank in the locus. The *Highest_ranked_gene_in_locus_plus_nearest* is similar, but prioritizes the nearest gene in case of ties. 


