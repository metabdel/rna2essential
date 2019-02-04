library(data.table)
library(readr)
library(tidyverse)


# load original data --- CCLE RNA-seq Linnorm
{
  # rna-seq data from DepQ2 release
  depQ4_cell <- fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/DepMap-2018q4-celllines.csv", 
                      data.table=FALSE, header=TRUE, stringsAsFactors = FALSE)
  rnaMapQ2 <- read_tsv("W:/Data/CCLE/broadInstitute/CCLE_DepMap_18Q2_RNAseq_HGNC_gene_map.txt") %>%
    dplyr::filter(locusType %in% c("gene with protein product", "protocadherin"))
  ccleLinnormQ2 <- fread("W:/Data/CCLE/broadInstitute/CCLE_DepMap_18Q2_RNAseq_count_20180502_linnorm_normalized.txt", 
                         data.table=FALSE, header=TRUE, stringsAsFactors = FALSE) %>%
    filter(ensembl %in% rnaMapQ2$originalSymbol) %>%
    arrange(ensembl) %>% 
    magrittr::set_rownames(rnaMapQ2$hgncSymbol) %>%
    select(-ensembl) %>%
    as.matrix()
  
  rnaQ2Cell <- filter(depQ4_cell, CCLE_Name %in% colnames(ccleLinnormQ2))
  ccleLinnormQ2 <- ccleLinnormQ2[, rnaQ2Cell$CCLE_Name] %>%
    magrittr::set_colnames(rnaQ2Cell$DepMap_ID)
  
  
  # rna-seq data from DepQ4 release (Linnorm normalized)
  rnaMapQ4 <- read_tsv("W:/Data/CCLE/broadInstitute/CCLE_depMap_18Q4_TPM_v2_gene_HGNC_map.txt") %>%
    dplyr::filter(locusType %in% c("gene with protein product", "protocadherin"))
  
  ccleLinnormQ4 <- readRDS("W:/Data/CCLE/broadInstitute/CCLE_depMap_18Q4_TPM_v2_Linnorm.rds") %>% 
    .[rnaMapQ4$original, ] %>%
    magrittr::set_rownames(rnaMapQ4$hgncSymbol) %>% 
    t()
  
  ccleRNAQ4 <- data.table::fread("W:/Data/CCLE/broadInstitute/CCLE_depMap_18Q4_TPM_v2.csv", data.table=FALSE, stringsAsFactors = FALSE) %>%
    tibble::column_to_rownames("V1") %>%
    as.matrix() %>%
    .[, rnaMapQ4$original] %>%
    magrittr::set_colnames(rnaMapQ4$hgncSymbol)
}

# load original data --- Avana 18Q4 data
{
  koMap <- fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/avana_depQ4_gene_dependency_gene_HGNC_map.txt",
                 data.table=FALSE, header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::filter(locusType %in% c("gene with protein product", "protocadherin"))
  
  geneDep <- fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/gene_dependency.csv", 
                   data.table=FALSE, header=TRUE, stringsAsFactors = FALSE) %>%
    tibble::column_to_rownames("DepMap_ID") %>%
    as.matrix() %>% 
    .[, koMap$oldName] %>%
    magrittr::set_colnames(koMap$hgncSymbol)
  
  geneEffect <- fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/gene_effect.csv",
                      data.table=FALSE, header=TRUE, stringsAsFactors = FALSE) %>%
    tibble::column_to_rownames("DepMap_ID") %>%
    as.matrix() %>%
    .[, koMap$oldName] %>%
    magrittr::set_colnames(koMap$hgncSymbol)
  
  commonCell_1 <- intersect(rownames(ccleLinnormQ4), rownames(geneDep))
  ccleLinnormQ4 <- ccleLinnormQ4[commonCell_1, ]
  geneDep <- geneDep[commonCell_1, ]
  geneEffect <- geneEffect[commonCell_1, ]
  
  repMap <- read_tsv("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/replicate_map_remove_batch3_control.txt") %>%
    setNames(c("replicate", "cell")) %>%
    dplyr::filter(cell %in% commonCell_1)
  guideMap <- read_tsv("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/guide_gene_HGNC_map_no_ambiguity.txt")
  guideMap <- guideMap[! duplicated(guideMap$sgrna), ]
  
  # ceres logFc data
  ceresFc <- fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/logfold_change.csv", 
                   data.table=FALSE, stringsAsFactors = FALSE) %>%
    tibble::column_to_rownames( "Construct Barcode") %>% 
    as.matrix() %>%
    .[, repMap$replicate]
  
}


# load processed data (selected genes for rnaseq-essentialGenes analysis) 
# include 
load("W:/scripts/rna2essential/DepQ4_select_gene_sgRNA_for_ccleLinnorm_prediction.RData")












