library(data.table)
library(readr)
library(tidyverse)
library(magrittr)


# expression/dependency data
{ 

  rnaMapQ2 <- readr::read_tsv("W:/projects/dependency_SL/ccleData/gene_annotation/CCLE_DepMap_18Q2_RNAseq_HGNC_gene_map.txt") %>%
    dplyr::filter(locusType %in% c("gene with protein product", "protocadherin")) %>%
    dplyr::arrange(originalSymbol)
  
  rnaQ2 <- readRDS("W:/Data/CCLE/broadInstitute/CCLE_DepMap_18Q2_RNAseq_reads_20180502_Linnorm.rds") %>%
    .[rnaMapQ2$originalSymbol, ] %>%
    magrittr::set_rownames(rnaMapQ2$hgncSymbol) %>%
    t()
  
  rnaMapQ4 <- read_tsv("W:/Data/CCLE/broadInstitute/CCLE_depMap_18Q4_TPM_v2_gene_HGNC_map.txt") %>%
    dplyr::filter(locusType %in% c("gene with protein product", "protocadherin"))
  
  rnaQ4 <- readRDS("W:/Data/CCLE/broadInstitute/CCLE_depMap_18Q4_TPM_v2.rds")  
  
  # read cn data
  cnMap <- readr::read_tsv("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/public_18Q4_gene_cn_gene_HGNC_map.txt")
  cn <- data.table::fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/public_18Q4_gene_cn.csv",
                          data.table=FALSE, stringsAsFactors = FALSE) %>%
    tibble::column_to_rownames("V1") %>%
    as.matrix() %>%
    .[, cnMap$original] %>%
    magrittr::set_colnames(cnMap$hgncSymbol) %>%
    apply(1, function(x) (2^x)*2) %>%
    t() %>%
    apply(2, function(x) ifelse(x < 0.25, 0, 1))  # if copy number less than 0.25, set the expression to 0
  
  cc <- intersect(rownames(cn), rownames(rnaQ4))
  cg <- intersect(colnames(cn), colnames(rnaQ4))
  cn <- cn[cc, cg]
  rnaQ4 <- rnaQ4[cc, cg] * cn
  
  # crispr knockout data
  koMap4 <- fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/avana_depQ4_gene_dependency_gene_HGNC_map.txt",
                 data.table=FALSE, header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::filter(locusType %in% c("gene with protein product", "protocadherin"))
  
  geneDep4 <- fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q4/gene_dependency.csv", 
                   data.table=FALSE, header=TRUE, stringsAsFactors = FALSE) %>%
    tibble::column_to_rownames("DepMap_ID") %>%
    as.matrix() %>% 
    .[, koMap$oldName] %>%
    magrittr::set_colnames(koMap$hgncSymbol)
  
  cellMapQ2 <- readr::read_tsv("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q2/18Q2_cell_map.txt")
  koMap2 <- fread("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q2/avana_18Q2_HGNC_gene_mapping.txt",
                  data.table=FALSE, header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::filter(locusType %in% c("gene with protein product", "protocadherin"))
  
  geneDepQ2 <- readRDS("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q2/gene_dependency.rds")
  geneEffect <- readRDS("W:/Data/genetic_screening/BroadInstitute/Broad_cancer_dependency_map_Avana_public_18Q2/gene_effect.rds")
  
  # validate cell line dependency and rna-seq
  geneDep5 <- geneDep4[! rownames(geneDep4) %in% cellMapQ2$Broad_ID, ]
  ccValidate <- intersect(rownames(geneDep5), rownames(rnaQ4))
  rnaQ4validate <- rnaQ4[ccValidate, ]
  geneDep5 <- geneDep5[ccValidate, ]
  
  save(geneDep5, rnaQ4validate, file="expression_dependency_for_validate.RData")
  # training cell line dependency and rna-seq
  cc2 <- Reduce(intersect, list(rownames(geneDepQ2), rownames(cn), rownames(geneDep4), rownames(rnaQ2)))
  cg2 <- intersect(colnames(cn), colnames(rnaQ4))
  cn <- cn[cc2, cg2]
  geneDep <- geneDep4[cc2, ]
  rnaQ4 <- rnaQ4[cc2, cg2]
  

  expSummary <- apply(rnaQ4, 2, function(x) 
    c(quantile(x, na.rm=TRUE), 
      variance=var(x, na.rm=TRUE),
      nonZero=sum(x >0, na.rm=TRUE))) %>% 
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    setNames(c("gene", "minimum", "quarter1","median", "quarter3", "maximum","variance", "nonZero"))
  
  expSummary2 <- filter(expSummary, nonZero >= 20)
  rnaQ4 <- rnaQ4[, expSummary2$gene]
  
  
  
  essentialGenes <- filter(combinedEssentiality, depValue > 0.75 | mageckPvalue < 0.05 | snMixture > 0.65)
  essentialCount <- table(essentialGenes$gene)
  essentialCount <- names(essentialCount[essentialCount > 150 & essentialCount < 300])
  geneDep2 <- gendDep[, essentialCount]
  
  depSummary <- apply(geneDep, 2, function(x) c(above_quarter3=sum(x > 0.75, na.rm = TRUE), below_quarter1=sum(x < 0.25, na.rm=TRUE))) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(new=above_quarter3 + below_quarter1)
  
  depSummary2 <- filter(depSummary, above_quarter3 > 50, below_quarter1 > 30, new > 250) %>%
    dplyr::filter(! gene %in% colnames(geneDep2))
  geneDep3 <- geneDep[, depSummary2$gene]
  
  save(rnaQ4, geneDep2,geneDep3, file="expression_TPM_4_gene_dependency.RData")
  
  
}

### select gene features

{
  ### pre-select gene variables
  depGenes3 <- colnames(geneDep3)
  sapply(depGenes3[1:60], function(x) {
    thisDep <- geneDep3[, x]
    preRNA <- vsurf_thres(rnaQ4, thisDep, x, eachSide=2000, useCores=30)
    saveRDS(preRNA, paste0(x, "_rna2d.rds"))
  })
  
  
  # select gene features by BART
  rna2dFiles <- list.files(pattern="*rna2d.rds")
  sapply(rna2dFiles, function(x) {
    gene <- gsub("_rna2d.rds", "", x)
    thisFc <- geneDep3[, gene]
    rna2d <- readRDS(x)
    # use bart_select3 or bart_select2 for selecting gene features
    start_time <- Sys.time()
    z <- bart_select2(rna2d, thisFc, useCores=28)
    end_time <- Sys.time()
    print(end_time - start_time)
    saveRDS(z, paste0("BART_selected_", gene, ".rds"))
    
  })
  
  
  # select gene features by bartMachine
  rna2dFiles <- list.files(pattern="*rna2d.rds")
  sapply(rna2dFiles, function(x) {
    gene <- gsub("_rna2d.rds", "", x)
    thisFc <- geneDep3[, gene]
    rna2d <- readRDS(x)
    z <- bartMachine_select2(rna2d, thisFc, useCores=28)
    saveRDS(z, paste0("bartMachine_selected_", gene, ".rds"))
    
  })
  
  
  # select gene features by xgboost
  rna2dFiles <- list.files(pattern="*rna2d.rds")
  sapply(rna2dFiles[c(5,11,17:104)], function(x) {
    gene <- gsub("_rna2d.rds", "", x)
    thisFc <- geneDep3[, gene]
    rna2d <- readRDS(x)
    z <- xgb_select(rna2d, thisFc, nthreads=25)
    saveRDS(z, paste0("xgboost_selected_", gene, ".rds"))
    
  })
  
  
  # select gene features by VSURF
  ### RRF select using rna2d selected by VSURF
  rna2dFiles <- list.files(pattern="*rna2d.rds")
  sapply(rna2dFiles, function(x) {
    gene <- gsub("_rna2d.rds", "", x)
    thisFc <- geneDep3[, gene]
    rna2d <- readRDS(x)
    z <- rf_VSURF_pred3(rna2d, thisFc, useCores=30)
    saveRDS(z, paste0("VSURF_selected_", gene, ".rds"))
    
  })
  
  # select gene features by ranger-mlr
  rna2dFiles <- list.files(pattern="*rna2d.rds")
  sapply(rna2dFiles, function(x) {
    gene <- gsub("_rna2d.rds", "", x)
    thisFc <- geneDep3[, gene]
    rna2d <- readRDS(x)
    z <- mlr_ranger_select(rna2d, thisFc, useCPUs=30)
    saveRDS(z, paste0("ranger_mlr_selected_", gene, ".rds"))
    
  })
  
  # select gene features by ranger
  rna2dFiles <- list.files(pattern="*rna2d.rds")
  sapply(rna2dFiles, function(x) {
    gene <- gsub("_rna2d.rds", "", x)
    thisFc <- geneDep3[, gene]
    rna2d <- readRDS(x)
    start_time <- Sys.time()
    z <- rf_ranger_select(rna2d, thisFc, useThreads=30)
    end_time <- Sys.time()
    print(end_time - start_time)
    saveRDS(z, paste0("ranger_selected_", gene, ".rds"))
    
  })
  
  
  
}







### train model using selected gene features
{
  load("expression_dependency_for_validate.RData")
  
  # bartMachine training for BART selected gene features
  BARTselects <- list.files(pattern="BART_selected*")
  sapply(BARTselects, function(x) {
    depGene <- gsub("BART_selected_", "", x)
    depGene <- gsub(".rds", "", depGene)
    rnaGene <- readRDS(x) %>% 
      lapply(function(y) y[1:15, ]) %>%
      purrr::reduce(full_join, by="gene") %>%
      tibble::column_to_rownames("gene") %>%
      apply(1, function(f) sum( ! is.na(f))) %>%
      sort(decreasing = TRUE) %>%
      .[1:15] %>%
      names()
    
    trainRNA <- rnaQ4[, rnaGene]
    trainFc <- geneDep3[, depGene]
    
    testRNA <- rnaQ4validate[, rnaGene]
    testFc <- geneDep5[, depGene]
    
    this_bartM_model <- bartMachine_train(trainRNA, trainFc, folds=10, useCPUs=30) 
    bartM_test <- bartMachine_predict(this_bartM_model$trainModel, testRNA, testFc)
    saveRDS(this_bartM_model, paste0("bartMachine_train_model_", depGene, ".rds"))
    saveRDS(bartM_test, paste0("bartMachine_test_", depGene, ".rds"))
    
  })
  
  # bartMachine-mlr training for BART selected gene features
  BARTselects <- list.files(pattern="BART_selected*")
  sapply(BARTselects, function(x) {
    depGene <- gsub("BART_selected_", "", x)
    depGene <- gsub(".rds", "", depGene)
    rnaGene <- readRDS(x) %>% 
      lapply(function(y) y[1:15, ]) %>%
      purrr::reduce(full_join, by="gene") %>%
      tibble::column_to_rownames("gene") %>%
      apply(1, function(f) sum( ! is.na(f))) %>%
      sort(decreasing = TRUE) %>%
      .[1:15] %>%
      names()
    
    trainRNA <- rnaQ4[, rnaGene]
    trainFc <- geneDep3[, depGene]
    
    testRNA <- rnaQ4validate[, rnaGene]
    testFc <- geneDep5[, depGene]
    
    this_bartM_model <- bartMachine_mlr_train(trainRNA, trainFc, useCPUs=30) 
    bartM_test <- bartMachine_mlr_predcit(this_bartM_model$trainModel, testRNA, testFc)
    saveRDS(this_bartM_model, paste0("bartMachine_mlr_train_model_", depGene, ".rds"))
    saveRDS(bartM_test, paste0("bartMachine_mlr_test_", depGene, ".rds"))
    
  })
  
  
  # bartMachine training for bartMachine selected gene features
  
  # bartMachine-mlr training for bartMachine selected gene features


  
  # xgboost train
  xgbSelects <- list.files(pattern="xgboost_selected*")
  sapply(xgbSelects, function(x) {
    depGene <- gsub("xgboost_selected_", "", x)
    depGene <- gsub(".rds", "", depGene)
    rnaGene <- readRDS(x) %>% 
      lapply(function(y)  
        lapply(y, function(z) z[1:15, c("Feature", "Gain")]) %>%
          purrr::reduce(full_join, by="Feature") )  %>%
      purrr:::reduce(full_join, by="Feature") %>%
      tibble::column_to_rownames("Feature") %>%
      apply(1, function(f) sum(! is.na(f))) %>%
      sort(decreasing = TRUE) %>%
      .[1:15] %>%
      names()
    
    trainRNA <- rnaQ4[, rnaGene]
    trainFc <- geneDep3[, depGene]
    
    testRNA <- rnaQ4validate[, rnaGene]
    testFc <- geneDep5[, depGene]
    
    xgbModel <- xgb_mlr_train(trainRNA, trainFc, useCPUs=20) 
    xgbTest <- xgb_mlr_predict(xgbModel$trainModel, testRNA, testFc)
    saveRDS(xgbModel, "xgboost_train_model_", depGene, ".rds")
    saveRDS(xgbTest, "xgboost_test_", depGene, ".rds")
    
  })
  
  
  
  # ranger-mlr train for ranger-mlr selected
  
  
  
  # ranger-mlr train for ranger selected
  rangerSelects <- list.files(pattern="ranger_selected*")
  sapply(rangerSelects, function(x) {
    depGene <- gsub("ranger_selected_", "", x)
    depGene <- gsub(".rds", "", depGene)
    rnaGene <- readRDS(x) %>% 
      lapply(function(y) y[1:15, ]) %>%
      purrr::reduce(full_join, by="gene") %>%
      tibble::column_to_rownames("gene") %>%
      apply(1, function(f) sum( ! is.na(f))) %>%
      sort(decreasing = TRUE) %>%
      .[1:15] %>%
      names()
    
    trainRNA <- rnaQ4[, rnaGene]
    trainFc <- geneDep3[, depGene]
    
    testRNA <- rnaQ4validate[, rnaGene]
    testFc <- geneDep5[, depGene]
    
    rangerModel <- ranger_mlr_train(trainRNA, trainFc, useCPUs=20) 
    rangerTest <- ranger_mlr_predict(rangerModel$trainModel, testRNA, testFc)
    saveRDS(xgbModel, "ranger_train_model_", depGene, ".rds")
    saveRDS(xgbTest, "ranger_test_", depGene, ".rds")
  })
    
}



### plot prediction and true 

