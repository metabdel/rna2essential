
# on cell server
# calculate the median value for each gene using selected sgrna and cell line replicates
fc_median <- function(foldChange, guideMap, repMap) {
  final <- sapply(unique(repMap$cell), function(x) {
    sapply(unique(guideMap$hgncSymbol), function(y) {
      thisRep <- filter(repMap, cell==x)
      thisSg <- filter(guideMap, hgncSymbol==y)
      thisFc <- foldChange[thisSg$sgrna, thisRep$replicate] %>%
        melt() %>%
        .$value %>%
        median(na.rm=TRUE)
    })
  })
}

# check sgRNA consistancy by median value 
remove_sgRNA <- function(p) {
  g <- unique(p$gene)
  leftSg <- lapply(g, function(x) {
    thisG <- filter(p, gene==x)
    mVec <- thisG$median
    names(mVec) <- thisG$sgrna
    mVec <- sort(mVec)
    minimum <- mVec[1]
    for(i in 1:(length(mVec)-1)) {
      if(mVec[length(mVec)]-minimum >= 1)
        mVec <- mVec[1:(length(mVec)-1)]
    }
    return(names(mVec))
  }) %>% unlist()
  
  leftSg <- filter(p, sgrna %in% leftSg)
  return(leftSg)
}


#

