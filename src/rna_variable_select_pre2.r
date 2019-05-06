library(caret)
library(magrittr)
library(mlr)
# initial variable selection using genes selected from ccleRNA by VSURF-threshold method

useGenes <- unique(selectGuide$hgncSymbol)
thisGene <- useGenes[1]
thisDep <- selectDep[thisGene, ]
thisFc <- selectGeneFc[thisGene, ]




# first step filtering of VSURF 
vsurf_thres <- function(rna, thisFc, testGene, eachSide=2000, useCores=28) {
  require(magrittr)
  rna <- rna[, ! colnames(rna) %in% testGene]
  # filter gene by correlation, choose the number of genes to keep at this stage
  thisCor <- apply(rna, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(eachSide)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(eachSide)
  selectCor <- c(posCor, minusCor)
  
  rna <- rna[, names(selectCor)]
  
  # use VSURF-threshold step to filter out genes
  start_time <- Sys.time()
  gene.thres <- VSURF::VSURF_thres(x=rna, y=thisFc, na.action=na.omit, parallel = TRUE, ncores=useCores)
  end_time <- Sys.time()
  print(end_time - start_time)
  # toy.interp <- VSURF::VSURF_interp(xInnerTrain, yInnerTrain, vars= toy.thres$varselect.thres, na.action=na.omit, parallel = TRUE, ncores=useCores)
  
  return(rna[, gene.thres$varselect.thres])
  
}


# boosting variable selection from xgboost package
xgb_select <- function(rna2d, thisFc, nthreads=28) {
  
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  
  outerList <- lapply(outerIndex, function(oi){
    xOuterTrain <- rna2d[oi, ]
    yOuterTrain <- thisFc[oi]
    
    
    innerIndex <- caret::createFolds(yOuterTrain, k=10, list=TRUE, returnTrain = TRUE)
    innerList <- lapply(innerIndex, function(ii){
      
      yInnerTrain <- yOuterTrain[ii]
      yInnerTest <- yOuterTrain[-ii]
      xInnerTrain <- xOuterTrain[ii, ]
      xInnerTest <- xOuterTrain[-ii, ]
      
      
      dtrain <- xgboost::xgb.DMatrix(data=xInnerTrain, label=yInnerTrain)
      dtest <- xgboost::xgb.DMatrix(data=xInnerTest, label=yInnerTest)
      watchlist <- list(train=dtrain, test=dtest)
      
      # eta <- 0
      # depth <- 0
      # rmse <- 100
      # for(thisDepth in seq(5, 20, 2)) {
      #  for(thisEta in seq(0.1, 0.9, 0.05)) {
      #    xgbSelect <- xgboost::xgb.train(data=dtrain, watchlist=watchlist, eta=thisEta, max.depth=thisDepth, nrounds=5000, early_stopping_rounds = 50, print_every_n = 500, nthread=nthreads, verbose=TRUE)
      #   thisRMSE <- xgbSelect$best_score
      #    if(thisRMSE < rmse) {
      #     rmse <- thisRMSE
      #      eta <- eta
      #      depth <- depth
      #    }
      #  }
      #}
      
      params <- list(booster="gbtree", objective="reg:linear", eta=0.01, gamma=0, max_depth=10)
      
      # xgbcv <- xgb.cv(params=params, data=dtrain, nrounds=1000, nfolds=10, stratified=TRUE, print_every_n=500, early_stopping_rounds=50, nthread=10, verbose=FALSE)
      xgbSelect <- xgboost::xgb.train(params=params, data=dtrain, watchlist=watchlist,
                                      nrounds=5000, early_stopping_rounds = 50, 
                                      print_every_n=1000, nthread=nthreads, verbose=0)
      varGain <- xgboost::xgb.importance(model=xgbSelect) %>%
        as.data.frame()
      
      # return(varGain)
    })
    
  })
  
  return(outerList)
  
}


# boosting variable selection from xgboost package
xgb_select2 <- function(rna2d, thisFc, nthreads=28) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  
  outerList <- lapply(outerIndex, function(oi){
    xOuterTrain <- rna2d[oi, ]
    xOuterTest <- rna2d[-oi, ]
    
    yOuterTrain <- thisFc[oi]
    yOuterTest <- thisFc[-oi]
      
    
    dtrain <- xgboost::xgb.DMatrix(data=xOuterTrain, label=yOuterTrain)
    dtest <- xgboost::xgb.DMatrix(data=xOuterTest, label=yOuterTest)
    watchlist <- list(train=dtrain, test=dtest)
      
    params <- list(booster="gbtree", objective="reg:linear", eta=0.01, gamma=0, max_depth=10)
      
    # xgbcv <- xgb.cv(params=params, data=dtrain, nrounds=1000, nfolds=10, stratified=TRUE, metrics=list("rmse"), print_every_n=500, early_stopping_rounds=50, nthread=10, verbose=FALSE)
    
    xgbSelect <- xgboost::xgb.train(params=params, data=dtrain, watchlist=watchlist,
                                    nrounds=5000, early_stopping_rounds = 50,
                                    print_every_n=1000, nthread=nthreads, verbose=0)
    varGain <- xgboost::xgb.importance(model=xgbSelect) %>%
      as.data.frame()
    
  })
  
  return(outerList)
  
}

# BART variable selection using variable selected from VSURF-first step 
bart_select <- function(rna2d, thisFc, useCores=28) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- rna2d[oi, ]
    yOuterTrain <- thisFc[oi]
    
    innerIndex <- caret::createFolds(yOuterTrain, k=10, list=TRUE, returnTrain = TRUE)
    
    # start_time <- Sys.time()
    innerList <- lapply(innerIndex, function(ii){
      yInnerTrain <- yOuterTrain[ii]
      yInnerTest <- yOuterTrain[-ii]
      xInnerTrain <- xOuterTrain[ii, ]
      xInnerTest <- xOuterTrain[-ii, ]
      
      dirichlet_select <- BART::mc.wbart(xInnerTrain, yInnerTrain, xInnerTest, cont=TRUE, a=0.5, rho=0.2*ncol(xInnerTrain), sigdf=3, sigquant=0.9, k=3, 
                                         sparse=TRUE, augment=TRUE, ntree=200, ndpost=8000, nskip=2000, printevery=1000, mc.cores = useCores, nice = 0)
      
      
      
      thisRMSE <- sqrt(mean((dirichlet_select$yhat.test.mean - yInnerTest)^2))
      
      theDF <- data.frame(gene=names(dirichlet_select$varprob.mean), varProb=dirichlet_select$varprob.mean, varCount=dirichlet_select$varcount.mean, stringsAsFactors = FALSE)
      
    })
    # end_time <- Sys.time()
    # about 7 minutes for inner loops
    
  })
  return(outerList)
  
}


bart_select2 <- function(rna2d, thisFc, useCores=28) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- rna2d[oi, ]
    xOuterTest <- rna2d[-oi, ]
    yOuterTrain <- thisFc[oi]
    yOuterTest <- thisFc[-oi]
    
    dirichlet_select <- BART::mc.wbart(xOuterTrain, yOuterTrain, xOuterTest, cont=TRUE, a=0.5, rho=0.2*ncol(xOuterTrain), sigdf=3, sigquant=0.9, k=3, 
                                       sparse=TRUE, augment=TRUE, ntree=200, ndpost=8000, nskip=2000, printevery=1000, mc.cores = useCores, nice = 0)
    
    # thisRMSE <- sqrt(mean((dirichlet_select$yhat.test.mean - yOuterTest)^2))
    
    theDF <- data.frame(gene=names(dirichlet_select$varprob.mean), varProb=dirichlet_select$varprob.mean, varCount=dirichlet_select$varcount.mean, stringsAsFactors = FALSE)
    
  })
  return(outerList)
  
}


bart_select3 <- function(rna2d, thisFc, useCores=28) {
  
  dirichlet_select <- BART::mc.wbart(rna2d, thisFc, cont=TRUE, a=0.5, rho=0.2*ncol(rna2d),
                                     sigdf=3, sigquant=0.9, k=3, sparse=TRUE, augment=TRUE, ntree=200,
                                     ndpost=8000, nskip=2000, printevery=1000, mc.cores = useCores, nice = 0)
  
  theDF <- data.frame(gene=names(dirichlet_select$varprob.mean), 
                      varProb=dirichlet_select$varprob.mean,
                      varCount=dirichlet_select$varcount.mean, 
                      stringsAsFactors = FALSE)
    
  return(theDF)
  
}


# bartMachine variable selection using variable selected from VSURF-first step 
bartMachine_select <- function(rna2d, thisFc, useCores=28) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- rna2d[oi, ]
    yOuterTrain <- thisFc[oi]
    
    # start_time <- Sys.time()
    bartMachine::set_bart_machine_num_cores(useCores)
    bart_machine <- bartMachine::bartMachine(as.data.frame(xOuterTrain), yOuterTrain)
    vs <- bartMachine::var_selection_by_permute_cv(bart_machine)
    vsSelect <- vs$important_vars_cv
    # end_time <- Sys.time()
    # print(end_time - start_time) #  about 10 minutes
    
  })
  return(outerList)
  
}



# bartMachine variable selection using variable selected from VSURF-first step 
bartMachine_select2 <- function(rna2d, thisFc, useCores=28) {
  # start_time <- Sys.time()
  bartMachine::set_bart_machine_num_cores(useCores)
  bart_machine <- bartMachine::bartMachine(as.data.frame(rna2d), thisFc)
  vs <- bartMachine::var_selection_by_permute_cv(bart_machine, k_folds=10)
  vsSelect <- vs$important_vars_cv
  # end_time <- Sys.time()
  # print(end_time - start_time) #  about 10 minutes
  
  return(vsSelect)
  
}


# ranger-random forest variable selection using variable selected from VSURF-first step 
rf_ranger_select <- function(rna2d, thisFc, useThreads=28) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- rna2d[oi, ]
    yOuterTrain <- thisFc[oi]
    
    # innerIndex <- caret::createFolds(ytrain, k=10, list=TRUE, returnTrain = TRUE)
    # start_time <- Sys.time()
    #innerMatrix <- sapply(innerIndex, function(ii){
    #  yInnerTrain <- ytrain[ii]
    #  yInnerTest <- ytrain[-ii]
    #  xInnerTrain <- xtrain[ii, ]
    #  xInnerTest <- xtrain[-ii, ]
      
      # random forest from ranger package
      
      rangerDF <- data.frame(xOuterTrain, response=yOuterTrain) 
      ranger_select <- ranger::ranger(formula=response~., data=rangerDF, num.trees=2000,
                                      max.depth=30, importance="impurity_corrected", 
                                      oob.error=TRUE, num.threads = useThreads)
      varImp <- ranger::importance(ranger_select) 
      varImp <- data.frame(gene=names(varImp), importance=varImp, stringsAsFactors = FALSE) %>%
        dplyr::arrange(desc(importance))
      # varImp <- varImp[order(names(varImp))]
  #  })
    
  })
  return(outerList)
  
}


# regularized random forest variable selection using RRF
rrf_select <- function(rna2d, thisFc) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- rna2d[oi, ]
    # xOuterTest <- rna2d[-oi, ]
    yOuterTrain <- thisFc[oi]
    # yOuterTest <- thisFc[-oi]
    innerIndex <- caret::createFolds(yOuterTrain, k=10, list=TRUE, returnTrain = TRUE)
    innerList <- parallel::mclapply(innerIndex, function(ii){
      yInnerTrain <- yOuterTrain[ii]
      yInnerTest <- yOuterTrain[-ii]
      xInnerTrain <- xOuterTrain[ii, ]
      xInnerTest <- xOuterTrain[-ii, ]
      
      rrf_fit <- RRF::RRF(x=xInnerTrain, y=yInnerTrain, xtest=xInnerTest, ytest=yInnerTest, ntree=2000, mtry=3, maxnodes=30, importance=TRUE, oob.prox=TRUE, flagReg=1)
      imp <- rrf_fit$importance
    }, mc.cores = 10)
  })
  return(outerList)
}

### time testing
start_time <- Sys.time()
rrfList <- rrf_select(rna2d, thisFc)
end_time <- Sys.time()

# random forest variable selection using caret-ranger
mlr_ranger_select <- function(rna2d, thisFc, useCPUs=20) {
  require(mlr)
  # library(parallelMap)
  parallelMap::parallelStartMulticore(useCPUs)
  
  # outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  # outerList <- lapply(outerIndex, function(oi){
    
  #  xOuterTrain <- rna2d[oi, ]
  #  yOuterTrain <- thisFc[oi]
    
  rangerDF <- data.frame(rna2d, response=thisFc) 
  selectTask <- makeRegrTask(data=rangerDF, target="response")
  ctrl <- makeFeatSelControlSequential(method="sfs", max.features = 30)
  learner <- makeLearner("regr.ranger", predict.type="response")
  inner <- makeResampleDesc(method="RepCV", reps=5, folds=10)
  selector <- selectFeatures(learner, selectTask, inner, measures = list(rmse, spearmanrho, mse), control=ctrl)
    
  # })
  
  parallelMap::parallelStop()
  # selectGenes <- sapply(outerList, function(y) y$x)
  return(selector)
  
}


# random forest variable selection using Boruta package
boruta_select <- function(rna2d, thisFc ) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- rna2d[oi, ]
    yOuterTrain <- thisFc[oi]
    
    boruta_fit <- Boruta::Boruta(x=data.frame(xOuterTrain), y=yOuterTrain, maxRuns=200, getImp=getImpRfZ)
    
  })
  return(outList)
}



# VSURF-random forest variable selection using VSURF package
rf_VSURF_pred <- function(rna2d, thisFc, testGene, useCores=28) {
  rna2d <- rna2d[, ! colnames(rna2d) %in% testGene]
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- rna2d[oi, ]
    yOuterTrain <- thisFc[oi]
    
    # library(VSURF)
    
    vsurf_step2 <- VSURF::VSURF_interp(x=xOuterTrain, y=yOuterTrain, vars=1:ncol(xOuterTrain), na.action = na.omit, parallel=TRUE, ncores=useCores)
    
    vsurf_step3 <- VSURF::VSURF_pred(x=xOuterTrain, y=yOuterTrain, err.interp =vsurf_step2$err.interp, vars=vsurf_step2$varselect.interp, na.action = na.omit, parallel=TRUE, ncores=useCores)
    
    trainVar <- colnames(xOuterTrain)
    varSelect <- list(varInterp=trainVar[vsurf_step2$varselect.interp], varPred=trainVar[vsurf_step3$varselect.pred])
    
  })
  return(outerList)
  
}

rf_VSURF_pred2 <- function(rna2d, thisFc, testGene, useCores=28) {
  rna2d <- rna2d[, ! colnames(rna2d) %in% testGene]
  # library(VSURF)
  # start_time <- Sys.time()
  vsurf_result <- VSURF::VSURF(x=rna2d, y=thisFc, na.action = na.omit, parallel=TRUE, ncores=useCores)
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  trainVar <- colnames(rna2d)
  varSelect <- list(varInterp=trainVar[vsurf_result$varselect.interp], varPred=trainVar[vsurf_result$varselect.pred])
  
  
  return(varSelect)
  
}

rf_VSURF_pred3 <- function(rna2d, thisFc, useCores=28) {
  
  start_time <- Sys.time()
  vsurf_step2 <- VSURF::VSURF_interp(x=rna2d, y=thisFc, vars=1:ncol(rna2d), na.action = na.omit, parallel=TRUE, ncores=useCores)
  
  vsurf_step3 <- VSURF::VSURF_pred(x=rna2d, y=thisFc, err.interp =vsurf_step2$err.interp, vars=vsurf_step2$varselect.interp, na.action = na.omit, parallel=TRUE, ncores=useCores)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  trainVar <- colnames(rna2d)
  varSelect <- list(varInterp=trainVar[vsurf_step2$varselect.interp], varPred=trainVar[vsurf_step3$varselect.pred])
  
  return(varSelect)
  
}








# caret xgboost variable selection using variable selected from VSURF-first step 
caret_xgb_select <- function(rna2d, thisFc, method="xgbTree", useCores=28) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- rna2d[oi, ]
    yOuterTrain <- thisFc[oi]
    
    
    library(doParallel)
    cl <- makeCluster(2)
    registerDoParallel(cl)
    subsets <- 5:50
    start_time <- Sys.time()
    ctrl <- caret::rfeControl(functions = caretFuncs, method = "cv", number = 10, verbose = FALSE, allowParallel = TRUE)
    rfe_profile <- caret::rfe(xOuterTrain, yOuterTrain, method="xgbTree", metric="RMSE",sizes=subsets, importance=TRUE, rfeControl = ctrl)
    end_time <- Sys.time()
    stopCluster(cl)
    
  })
  
  return(outerList)
  
}



# test lasso
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


### parameter grid search of BART
innerGrid <- function(xtrain, ytrain, useCores=15) {
  
  rho <- c(0.2, 0,5, 0.8)
  a <- c(0.5, 0.7)
  sigdf <- c(2, 3)
  sigquant <- c(0.9, 0.95)
  k <- c(2, 3, 4, 5)
  
  
  innerIndex <- caret::createFolds(ytrain, k=10, list=TRUE, returnTrain = TRUE)
  
  innerList <- lapply(innerIndex, function(ii){ 
    
    yInnerTrain <- ytrain[ii]
    yInnerTest <- ytrain[-ii]
    xInnerTrain <- xtrain[ii, ]
    xInnerTest <- xtrain[-ii, ]
    # innerDepTrain <- depTrain[ii]
    # innerDepTest <- depTrain[-ii]
    
    thisCor <- apply(xInnerTrain, 2, function(y) cor(y, yInnerTrain, method="spearman"))
    thisCor <- thisCor[! is.na(thisCor)]
    
    posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
    minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
    selectCor <- c(posCor, minusCor)
    
    xInnerTrain <- xInnerTrain[, names(selectCor)]
    thisVariables <- c()
    thisRMSE <- 100
    thisParameters <- c()
    
    for(l1 in a) {
      for(l2 in rho) {
        for(l3 in sigdf) {
          for(l4 in sigquant) {
            for(l5 in k) {
              dirichlet_select <- mc.wbart(xInnerTrain, yInnerTrain, xInnerTest, cont=TRUE, a=l1, rho=l2*1000, sigdf=l3, sigquant=l4, k=l5, 
                                           sparse=TRUE, augment=TRUE, ntree=200, ndpost=10000, nskip=5000, printevery=1000, mc.cores = useCores)
              
              rmse <- sqrt(mean((dirichlet_select$yhat.test.mean - yInnerTest)^2))
              
              if (rmse <= thisRMSE) {
                thisRMSE <- rmse
                thisVariables <- dirichlet_select$varprob.mean
                thisParameters <- c(a=l1, rho=l2*1000, sigdf=l3, sigquant=l4, k=l5)
              }
              
            }
          }
        }
      }
    }
    
    return(list(thisVariables, thisRMSE, thisParameters))
    
    
  })
  
  return(innerList)
  
}

# inner cross-validation of blasso
innerBlasso <- function(xtrain, ytrain) {
  innerIndex <- caret::createFolds(ytrain, k=10, list=TRUE, returnTrain = TRUE)
  
  #thisVariables <- c()
  #rmse <- c()
  
  varList <- lapply(innerIndex, function(ii){
    #for(ii in innerIndex){ 
    
    yInnerTrain <- ytrain[ii]
    yInnerTest <- ytrain[-ii]
    xInnerTrain <- xtrain[ii, ]
    xInnerTest <- xtrain[-ii, ]
    # innerDepTrain <- depTrain[ii]
    # innerDepTest <- depTrain[-ii]
    
    thisCor <- apply(xInnerTrain, 2, function(y) cor(y, yInnerTrain, method="spearman"))
    thisCor <- thisCor[! is.na(thisCor)]
    
    posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
    minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
    selectCor <- c(posCor, minusCor)
    
    xInnerTrain <- xInnerTrain[, names(selectCor)]
    
    
    
    library(doParallel)
    cl <- makeCluster(8)
    registerDoParallel(cl)
    a <- seq(0.1, 0.9, 0.05)
    search <- foreach(i=a, .combine=rbind) %dopar% {
      cv <- glmnet::cv.glmnet(xInnerTrain, yInnerTrain, family="gaussian", nfold=10, alpha=i, type.measure = "deviance", parallel = TRUE)
      data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
    }
    cv3 <- search[search$cvm == min(search$cvm), ]
    stopCluster(cl)
    
    elastic_fit <- glmnet(xInnerTrain, yInnerTrain, family="gaussian", lambda = cv3$lambda.1se, alpha = cv3$alpha)
    lasso.coef <- predict(elastic_fit, type="coefficients", s=cv$lambda.lse)[,1]
    lasso.coef <- lasso.coef[lasso.coef != 0]
    lasso.coef <- lasso.coef[2:length(lasso.coef)]
    
    
  })
  return(varList)
}

# caret recursive feature elimination example
{
  library(caret)
  
  data(mdrr)
  mdrrDescr <- mdrrDescr[,-nearZeroVar(mdrrDescr)]
  mdrrDescr <- mdrrDescr[, -findCorrelation(cor(mdrrDescr), .8)]
  
  set.seed(1)
  inTrain <- createDataPartition(mdrrClass, p = .75, list = FALSE)[,1]
  
  train <- mdrrDescr[ inTrain, ]
  test  <- mdrrDescr[-inTrain, ]
  trainClass <- mdrrClass[ inTrain]
  testClass  <- mdrrClass[-inTrain]
  
  set.seed(2)
  ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 5, verbose = FALSE)
  
  rf_profile <- rfe(train, trainClass,
                    ntree = 50,
                    rfeControl = ctrl)
  
}








# train and test should be sampled based on the distribution of response (y) value

train <- sample(1:nrow(ccleRNA), round(nrow(ccleRNA)*0.8))
x <- ccleRNA[, names(selectCor)]
xtrain <- x[train, ]
xtest <- x[-train, ]
ytrain <- thisFc[train]
ytest <- thisFc[-train]
### BART variable selection

# GSE method to perform variable selection
gseSelect <- mc.wbart.gse(xtrain, ytrain, ntree=30, P=100, mc.cores=15)
vsGenes <- gseSelect$which
vsGenes <- selectCor[vsGenes]



