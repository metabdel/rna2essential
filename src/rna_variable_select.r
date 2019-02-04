# for each gene, set the training and testing samples

# LASSO PART
selectGenes <- colnames(selectBCMI)
# parallel::mclapply(selectGenes, function(y) {

# prepare train and test rna data 
depValue <- geneDep2[, selectGenes[1]]

pos <- depValue[depValue > 0.75]
posSample <- sample(1:length(pos), round(length(pos)*0.75), replace=FALSE)
posTrain <- pos[posSample]
posTest <- pos[ -posSample]

neg <- depValue[depValue < 0.25]
negSample <- sample(1:length(neg), round(length(neg)*0.75), replace=FALSE)
negTrain <- neg[negSample]
negTest <- neg[ -negSample]

rest <- depValue[depValue >= 0.25 & depValue <= 0.75]
if(length(rest) > 100) { 
  rest <- rest[sample(1:length(rest),100, replace=FALSE)]
  restSample <- sample(1:length(rest), 50, replace=FALSE)
  restTrain <- rest[restSample]
  restTest <- rest[-restSample]
} else {
  restSample <- sample(1:length(rest), round(length(rest)*0.75), replace=FALSE)
  restTrain <- rest[restSample]
  restTest <- rest[-restSample]
}

train <- c(posTrain, negTrain, restTrain)
test <- c(posTest, negTest, restTest)
depTrain <- depValue[names(train)]
depTest <- depValue[names(test)]

rnaTrain <- ccleLinnorm[names(train), ]
rnaTest <- ccleLinnorm[names(test), ]

# rnaGenes <- rownames(selectBCMI)
#rnaGenes <- colnames(ccleLinnorm)
rnaGenes <- rnaSelect$gene
thisPvalue <- sapply(rnaGenes, function(x) {
  
  
  thisTrain <- rnaTrain[, x] %>% 
    as.data.frame() %>%
    setNames(x) %>%
    dplyr::mutate(depValue=depTrain)
  
  thisTest <- rnaTest[, x] %>%
    as.data.frame() %>%
    setNames(x) %>%
    dplyr::mutate(depValue=depTest)
  #testActual <- ifelse(thisTest$depValue >= 0.65, "e", "ne")
  
  glm.fits <- glm(depValue~., data=thisTrain, family=binomial)
  pValue <- summary(glm.fits)$coef
  pValue <- ifelse(nrow(pValue) ==2, pValue[2,4], 1)
  #  glm.probs <- predict(glm.fits, thisTest, type="response")
  #  glm.pred <- rep("ne", nrow(thisTest))
  #  glm.pred[glm.probs > 0.75] <- "e"
  #  return(mean(glm.pred==testActual))
  #} else{
  #  return(0)
  #}
  
})


thisPvalue <- names(sort(thisPvalue[thisPvalue < 0.05]))




thisTrain <- rnaTrain[, rownames(rnaTrain) %in% geneSet] %>% 
  as.data.frame() %>%
  setNames(geneSet) %>%
  dplyr::mutate(depValue=depTrain)

thisTest <- rnaTest[, rownames(rnaTest) %in% geneSet] %>%
  as.data.frame() %>%
  setNames(geneSet) %>%
  dplyr::mutate(depValue=depTest)

glm.fits <- glm(depValue~., data=thisTrain, family=binomial)
summary(glm.fits)









library(glmnet)
grid =10^ seq (10,-2, length =100)
selectGenes <- colnames(selectBCMI)


thisGene <- selectGenes[1]
yfc <- selectGeneFc[, thisGene]
miValue <- sort(selectBCMI[, thisGene], decreasing = TRUE)

top10 <- names(z)

rnaSample <- ccleLinnorm[, top10]
train <- sample(1:nrow(ccleLinnorm), round(nrow(ccleLinnorm)*0.75))
lasso.mod <- glmnet(rnaSample[train,], yfc[train], alpha=1, lambda=grid)

require(doMC)
registerDoMC(cores=2)


cv.out <- cv.glmnet(rnaSample[train,], yfc[train], alpha=1)
bestlam =cv.out$lambda.min
lasso.pred <- predict(lasso.mod, s=cv.out$lambda.min, newx=rnaSample[-train,])
mean((lasso.pred[,1]-yfc[-train])^2)

out <- glmnet(rnaSample, yfc, alpha=1, lambda=grid)
lasso.coef <- predict(out, type="coefficients", s=bestlam)[1:32,]

train <- sample(1:nrow(ccleLinnorm), round(nrow(ccleLinnorm)*0.75), replace = FALSE)
#rnaSample <- ccleLinnorm[, start:end] 

testLasso <- function(ccleLinnorm, train, yfc, interval=10, alpha=1) {
  iterationNo <- ceiling(ncol(ccleLinnorm)/interval)
  out <- lapply(1:iterationNo, function(idx) {
    
    out2 <- tryCatch (
      {
        start <- (idx-1)*interval + 1
        end <- min(idx*interval, 17976)
        grid =10^ seq (10,-2, length =100)
        
        #train <- sample(1:nrow(ccleLinnorm), round(nrow(ccleLinnorm)*0.75))
        rnaSample <- ccleLinnorm[, start:end] 
        
        cv.out <- cv.glmnet(rnaSample[train,], yfc[train], alpha=alpha)
        bestlam =cv.out$lambda.min
        out <- glmnet(rnaSample, yfc, alpha=1, lambda=grid)
        lasso.coef <- predict(out, type="coefficients", s=bestlam)[2:interval,]
        c <- lasso.coef[lasso.coef != 0]
        return(c)
      },
      error=function(e) { return(NULL)}
    )
    return(out2)
  }) %>% unlist()
  
}

# for one dependency gene
c0 <- testLasso(ccleLinnorm[, names(y)], train, yfc, 25, 1) 
c1 <- testLasso(ccleLinnorm, train, yfc, 10, 1) 

# for multiple dependency genes
c <- sapply(colnames(ccleLinnorm), function(x) cor(ccleLinnorm[, x], yfc, method="spearman"))


# 
train <- sample(1:nrow(ccleLinnorm), round(nrow(ccleLinnorm)*0.8))





# BART part
sparseSelect <- mc.wbart(xtrain, ytrain, sparse=TRUE, augment = TRUE, cont=TRUE, ntree=100, ndpost=5000, nskip=3000, printevery = 1000, transposed = TRUE, mc.cores=15)
gseSelect <- mc.wbart.gse(xtrain, ytrain, ntree=30, P=100, ndpost=8000, nskip=3000, printevery = 1000,  mc.cores=15)

vsGenes <- gseSelect$which
bfp1 <- wbart(xtrain[,s], ytrain,  sparse=TRUE, ntree=200, ndpost=6000, nskip=1000)
pred <- predict(bfp1, xtest[, s])


# BART variable selection
{
  useGenes <- unique(selectGuide$hgncSymbol)
  #sigmaList <- list()
  vsGeneList <- list()
  yhatTestList <- list()
  yhatTestMeanList <- list()
  selectCorList <- list()
  train <- sample(1:nrow(ccleLinnorm), round(nrow(ccleLinnorm)*0.8))
  for (i in 1:length(useGenes)) {
    
    thisGene <- useGenes[i]
    thisFc <- selectGeneFc[thisGene, ]
    
    thisCor <- apply(ccleLinnorm, 2, function(x) cor(x, thisFc, method="spearman"))
    thisCor <- thisCor[! is.na(thisCor)]
    
    posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(300)
    minusCor <- thisCor[thisCor <0] %>% sort() %>% head(300)
    selectCor <- c(posCor, minusCor)
    #selectCorList[[i]] <- selectCor[vsGeneList[[i]]] 
  
    x <- ccleLinnorm[, names(selectCor)]
    xtrain <- x[train, ]
    xtest <- x[-train, ]
    ytrain <- thisFc[train]
    ytest <- thisFc[-train]
    
    # bfSelect <- mc.wbart(xtrain, ytrain, sparse=TRUE, augment = TRUE, ntree=50, ndpost=5000, nskip=2000, mc.cores=10)
    testGSE <- mc.wbart.gse(xtrain, ytrain, ntree=30, P=100, mc.cores=6)
    vsGenes <- testGSE$which
    bfp <- mc.wbart(xtrain[,vsGenes], ytrain, x.test = xtest[, vsGenes], ntree=200, ndpost=6000, nskip=2000, mc.cores = 10)
    
    vsGeneList[[i]] <- selectCor[vsGenes]
    #sigmaList[[i]] <- bfSelect$sigma
    yhatTestList[[i]] <- bfp$varprob
    yhatTestMeanList[[i]] <- bfp$varprob.mean
    
  }
  
}



# BART test

train <- sample(1:nrow(ccleRNA), round(nrow(ccleRNA)*0.8))

bartSelect <- function(ccleRNA, guideGenes, selectGeneFc, useCores=10, selectNo=300) {
  
  selectList <- lapply(guideGenes, function(i) {
    
    yfc <- selectGeneFc[thisGene, ]
    thisCor <- apply(ccleRNA, 2, function(x) cor(x, yfc, method="spearman"))
    thisCor <- thisCor[! is.na(thisCor)]
    
    posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(selectNo)
    minusCor <- thisCor[thisCor < 0] %>% sort() %>% head(selectNo)
    selectCor <- c(posCor, minusCor)
    
    if(length(selectCor) > 0 ) {
      x <- ccleLinnorm[, names(selectCor)]
      xtrain <- x[train, ]
      ytrain <- yfc[train]
      testGSE <- mc.wbart.gse(xtrain, ytrain, ntree=30, P=100, mc.cores=useCores)
      z <- testGSE$which
      z <- selectCor[z]
      return(z)
    } else {
      return(NULL)
    }
    
  })
  return(selectList)
}


s <- c("RHOB", "TSPAN12", "CACNA1H", "B3GAT1", "NKD2", "PBX1", "SELENBP1", "GPR37", "SYT8", "IL12RB1",  "XAF1", "CRACR2A", "F8", 
       "CCNA1", "EVI2B", "MNX1", "REEP1", "HOOK1","BRSK2", "HES4", "MYH14", "CNPY3", "CGAS", "SCIN", "AXIN2", "FAM47E", "RHOU",
       "MLXIPL",  "PDZRN3", "CGREF1", "CARD16", "C10orf55", "EVI2A", "ST3GAL6",  "ANKLE1")




### bartMachine part

bart_machine <- bartMachine(as.data.frame(xtrain), ytrain)
vs <- var_selection_by_permute(bart_machine, bottom_margin=10, num_permute_samples=10)
vs <- var_selection_by_permute_cv(bart_machine)

vsSelect <- vs$important_vars_cv
bart_machine2 <- bartMachine(as.data.frame(xtrain[, vsSelect]), ytrain)
vs2 <- var_selection_by_permute_cv(bart_machine)

bart_machine_cv <- bartMachineCV(as.data.frame(xtrain[, vsSelect]), ytrain)

oos_perf <- bart_predict_for_test_data(bart_machine_cv, as.data.frame(xtest[, vs$important_vars_cv]), ytest)
oos_perf2 <- bart_predict_for_test_data(bart_machine_cv2, as.data.frame(xtest[, vs2$important_vars_cv]), ytest)

# bartMachine variable selection
{
  useGenes <- unique(selectGuide$hgncSymbol)
  vsList <- list()
  oos_perf_List <- list()
  
  for (i in 29:length(useGenes)) {
    
    thisGene <- useGenes[i]
    
    thisFc <- selectGeneFc[thisGene, ]
    #thisDep <- selectDep[thisGene, ]
    thisCor <- apply(ccleLinnorm, 2, function(x) cor(x, thisFc, method="spearman"))
    thisCor <- thisCor[! is.na(thisCor)]
    
    posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(300)
    minusCor <- thisCor[thisCor <0] %>% sort() %>% head(300)
    selectCor <- c(posCor, minusCor)
    
    
    x <- ccleLinnorm[, names(selectCor)]
    bart_machine <- bartMachine(as.data.frame(x), thisFc)
    vs <- var_selection_by_permute_cv(bart_machine)
    vsSelect <- vs$important_vars_cv
    if(length(vsSelect) == 0) {
      vsList[[i]] <- NULL
      oos_perf_List[[i]] <- NULL
    } else {
      train <- sample(1:nrow(ccleLinnorm), round(nrow(ccleLinnorm)*0.8))
      xtrain <- x[train,]
      ytrain <- thisFc[train]
      xtest <- x[-train, ]
      ytest <- thisFc[-train]
      bart_machine_cv <- bartMachineCV(as.data.frame(xtrain[, vsSelect]), ytrain)
      
      oos_perf <- bart_predict_for_test_data(bart_machine_cv, as.data.frame(xtest[, vsSelect]), ytest)
      
      vsList[[i]] <- vsSelect
      oos_perf_List[[i]] <- oos_perf
    }
  }
  
}




### chaotic test 
thisGene <- "FERMT2"
thisFc <- selectGeneFc[thisGene, ]
thisDep <- selectDep[thisGene, ]
thisCor <- apply(ccleLinnorm, 2, function(x) cor(x, thisFc, method="spearman"))
thisCor <- thisCor[! is.na(thisCor)]

posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
selectCor <- c(posCor, minusCor)


x <- ccleLinnorm[, names(selectCor)]
xtrain <- x[train, ]
xtest <- x[-train, ]
ytrain <- thisFc[train]
ytest <- thisFc[-train]








thisSg <- filter(guideEfficacy, hgncSymbol=="PTEN") %>% .$sgrna
thisRep <- filter(hartCell, cell=="ACH-000811") %>% .$replicate
thisFc <- chooseFc[thisSg, ] %>% t()
qplot(thisFc, geom="density")

cor(thisFc, method="spearman")
boxplot(thisFc, las=2)
thisFc2 <- ceresFc[thisSg, thisRep] %>% t()
thisFc2

# test1 <- bartSelect(ccleLinnorm, g[1:3], selectGeneFc)
# test2 <- bartSelect(ccleRNA, g[1:3], selectGeneFc)



# hart-essential genes that actually are not essential in most cell lines
c("CDK17", "MSANTD3", "WDR60", "THOP1", "ZBTB48", "C12orf66", "EIF5B")

# 





testForest <- function(ccleLinnorm, idx) {
  start <- (idx-1)*10 + 1
  end <- min(idx*10, 17976)
  
  
}

testSVM <- function(ccleLinnorm, yfc, idx) {
  
}


c <- parallel::mclapply(1:ncol(ccleLinnorm), function(x) {
  out <- tryCatch(
    {mi <- mpmi::cmi.pw(ccleLinnorm[,x], yfc)$bcmi
    return(mi)
    },
    error=function(e) return(NULL)
  )
  return(out)
}, mc.cores=10)


