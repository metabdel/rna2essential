library(BART)
library(magrittr)
library(caret)
library(magrittr)


useGenes <- unique(selectGuide$hgncSymbol)
thisGene <- useGenes[1]
thisDep <- selectDep[thisGene, ]
thisFc <- selectGeneFc[thisGene, ]

outerIndex <- createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)

ytrain <- thisFc[oi]
ytest <- thisFc[-oi]
xtrain <- ccleRNA[oi, ]
xtest <- ccleRNA[-oi, ]

# BART variable selection
bart_select <- function(ccleRNA, thisFc, useCores=20) {
  
  thisCor <- apply(ccleRNA, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
  selectCor <- c(posCor, minusCor)
  
  ccleRNA <- ccleRNA[, names(selectCor)]
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    yOuterTrain <- thisFc[oi]
    
    innerIndex <- caret::createFolds(yOuterTrain, k=10, list=TRUE, returnTrain = TRUE)
    
    # start_time <- Sys.time()
    innerList <- lapply(innerIndex, function(ii){
      yInnerTrain <- yOuterTrain[ii]
      yInnerTest <- yOuterTrain[-ii]
      xInnerTrain <- xOuterTrain[ii, ]
      xInnerTest <- xOuterTrain[-ii, ]
      
      dirichlet_select <- BART::mc.wbart(xInnerTrain, yInnerTrain, cont=TRUE, a=0.5, rho=0.2*1000, sigdf=3, sigquant=0.9, k=3, 
                                         sparse=TRUE, augment=TRUE, ntree=200, ndpost=10000, nskip=2000, printevery=1000, mc.cores = useCores, nice = 0)
      
      thisRMSE <- sqrt(mean((dirichlet_select$yhat.test.mean - yInnerTest)^2))
      
      return(list(varProb=dirichlet_select$varprob.mean, varCount=dirichlet_select$varcount.mean, rmse=thisRMSE))
      
    })
    # end_time <- Sys.time()
    
    # return(innerList)
  })
  return(outerList)
  
}



bartResult <- bart_select(ccleRNA, thisFc, useCores=20)
bartResult <- lapply(bartResult, function(x) 
  lapply(x, function(x) sort(x$varCount, decreasing = TRUE) %>% 
           head(50) %>% 
           names()) %>%
    unlist() ) %>% 
  unlist() %>%
  talbe() %>%
  sort(decreasing = TRUE) 




# BART predict
bart_predict <- function(x, y, bart_select_var) {
  
}


# bartMachine variable selection
bartMachine_select <- function(xtrain, ytrain, useCores=30) {
  
  thisCor <- apply(ccleRNA, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
  selectCor <- c(posCor, minusCor)
  
  ccleRNA <- ccleRNA[, names(selectCor)]
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    yOuterTrain <- thisFc[oi]
    
    innerIndex <- caret::createFolds(yOuterTrain, k=10, list=TRUE, returnTrain = TRUE)
    
    # start_time <- Sys.time()
    innerList <- lapply(innerIndex, function(ii){
      yInnerTrain <- yOuterTrain[ii]
      yInnerTest <- yOuterTrain[-ii]
      xInnerTrain <- xOuterTrain[ii, ]
      xInnerTest <- xOuterTrain[-ii, ]
      
      set_bart_machine_num_cores(useCores)
      bart_machine <- bartMachine::bartMachine(as.data.frame(xInnerTrain), yInnerTrain)
      vs <- bartMachine::var_selection_by_permute_cv(bart_machine)
      vsSelect <- vs$important_vars_cv
      
    })
    #end_time <- Sys.time()
    # return(innerList)
  })
  return(outerList)
  
}


# random forest variable selection using ranger package
rf_ranger_select <- function(ccleRNA, thisFc, useThreads=30) {
  
  thisCor <- apply(ccleRNA, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
  selectCor <- c(posCor, minusCor)
  
  ccleRNA <- ccleRNA[, names(selectCor)]
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    yOuterTrain <- thisFc[oi]
    
    innerIndex <- caret::createFolds(ytrain, k=10, list=TRUE, returnTrain = TRUE)
    # start_time <- Sys.time()
    innerMatrix <- sapply(innerIndex, function(ii){
      yInnerTrain <- ytrain[ii]
      yInnerTest <- ytrain[-ii]
      xInnerTrain <- xtrain[ii, ]
      xInnerTest <- xtrain[-ii, ]
      
      # random forest from ranger package
      rangerDF <- data.frame(xInnerTrain, response=yInnerTrain) 
      ranger_select <- ranger::ranger(formula=response~., data=rangerDF, num.trees=2000, max.depth=60, importance="impurity_corrected", oob.error=TRUE, num.threads = useThreads)
      varImp <- ranger::importance(ranger_select)
      varImp <- varImp[order(names(varImp))]
    })
    #end_time <- Sys.time()
    
    # return(innerMatrix)
  })
  return(outerList)
  
}


# random forest variable selection using VSURF package
rf_VSURF_select <- function(ccleRNA, thisFc, useCores=26) {
  
  thisCor <- apply(ccleRNA, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
  selectCor <- c(posCor, minusCor)
  
  ccleRNA <- ccleRNA[, names(selectCor)]
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    yOuterTrain <- thisFc[oi]
    
    innerIndex <- caret::createFolds(ytrain, k=10, list=TRUE, returnTrain = TRUE)
    
    innerMatrix <- sapply(innerIndex, function(ii){
      yInnerTrain <- ytrain[ii]
      yInnerTest <- ytrain[-ii]
      xInnerTrain <- xtrain[ii, ]
      xInnerTest <- xtrain[-ii, ]
      
      library(VSURF)
      vsurfSelect <- VSURF::VSURF(x=xInnerTrain, y=yInnerTrain, na.action = na.omit, parallel=TRUE, ncores=useCores)
      trainVar <- colnames(xInnerTrain)
      varSelect <- trainVar[vsurfSelect$varselect.pred]
      # toy.thres <- VSURF::VSURF_thres(xInnerTrain, yInnerTrain, na.action=na.omit, parallel = TRUE, ncores=useCores)
      # toy.interp <- VSURF::VSURF_interp(xInnerTrain, yInnerTrain, vars= toy.thres$varselect.thres, na.action=na.omit, parallel = TRUE, ncores=useCores)
      
    })
    # return(innerMatrix)
  })
  return(outerList)
  
}





# boosting variable selection from gbm package
gbm_select <- function(ccleRNA, thisFc) {
  
  # filter ccleRNA by correlation
  thisCor <- apply(ccleRNA, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
  selectCor <- c(posCor, minusCor)
  
  ccleRNA <- ccleRNA[, names(selectCor)]
  
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    yOuterTrain <- thisFc[oi]
    
    
    innerIndex <- caret::createFolds(yOuterTrain, k=10, list=TRUE, returnTrain = TRUE)
    innerMatrix <- sapply(innerIndex, function(ii){
      
      yInnerTrain <- yOuterTrain[ii]
      yInnerTest <- yOuterTrain[-ii]
      xInnerTrain <- xOuterTrain[ii, ]
      xInnerTest <- xOuterTrain[-ii, ]
      
      # boostDF <- data.frame(xInnerTrain, response=yInnerTrain)
      # gbmSelect <- gbm::gbm(formula=response~., distribution = "gaussian", data=boostDF, n.trees=1000, cv.folds = 10, n.cores=useCores)
      gbmSelect <- gbm::gbm.fit(x=xInnerTrain, y=yInnerTrain, distribution = "gaussian", n.trees=5000, shrinkage = 0.01, verbose=FALSE)
      
      gbmSummary <- summary(gbmSelect, plotit=FALSE, mormalize=FALSE, order=FALSE)
      varInfluence <- gbmSummary$rel.inf 
      names(varInfluence) <- gbmSummary$var
      return(varInfluence)
    })
    # return(innerMatrix)
  })
  
  return(outerList)
  
}


# boosting variable selection from xgboost package
xgb_select <- function(ccleRNA, thisFc, nthreads=20) {
  
  # filter ccleRNA by correlation
  thisCor <- apply(ccleRNA, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
  selectCor <- c(posCor, minusCor)
  
  ccleRNA <- ccleRNA[, names(selectCor)]
  
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
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
      
      eta <- 0
      depth <- 0
      rmse <- 100
      for(eta in seq(5, 20, 2)) {
        for(depth in seq(0.1, 0.9, 0.05)) {
          xgbSelect <- xgboost::xgb.train(data=dtrain, watchlist=watchlist, eta=eta, max.depth=depth, nrounds=5000, early_stopping_rounds = 50, print_every_n = 500, nthread=nthreads, verbose=TRUE)
          if(xgbSelect$best_score < rmse) {
            rmse <- xgbSelect$best_score
            eta <- eta
            depth <- depth
          }
        }
      }
      
      xgbSelect <- xgboost::xgb.train(data=dtrain, watchlist=watchlist, eta=eta, max.depth=depth, nrounds=5000, early_stopping_rounds = 50, print_every_n = 500, nthread=nthreads, verbose=TRUE)
      varGain <- xgboost::xgb.importance(model=xgbSelect) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("Feature") %>%
        as.matrix() %>%
        .[, 1]
      
      # return(varGain)
    })
    # return(innerList)
  })
  
  return(outerList)
  
}

# random forest variable selection using party package
rfSelect_3 <- function(ccleRNA, thisFc, useThreads=20) {
  
  thisCor <- apply(ccleRNA, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
  selectCor <- c(posCor, minusCor)
  
  ccleRNA <- ccleRNA[, names(selectCor)]
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    yOuterTrain <- thisFc[oi]
    
    innerIndex <- caret::createFolds(ytrain, k=10, list=TRUE, returnTrain = TRUE)
    
    innerMatrix <- sapply(innerIndex, function(ii){
      yInnerTrain <- ytrain[ii]
      yInnerTest <- ytrain[-ii]
      xInnerTrain <- xtrain[ii, ]
      xInnerTest <- xtrain[-ii, ]
      
      
      library(party)
      partyDF <- data.frame(xInnerTrain, response=yInnerTrain)
      partySelect <- party::cforest(response ~., data=partyDF, control=cforest_unbiased(ntree=1000))
      
      
    })
    # return(innerMatrix)
  })
  return(outerList)
  
}


# elastic net variable selection
elasticSelect <- function(ccleRNA, thisFc) {
  
  thisCor <- apply(ccleRNA, 2, function(y) cor(y, thisFc, method="spearman"))
  thisCor <- thisCor[! is.na(thisCor)]
  
  posCor <- thisCor[thisCor > 0] %>% sort(decreasing = TRUE) %>% head(500)
  minusCor <- thisCor[thisCor <0] %>% sort() %>% head(500)
  selectCor <- c(posCor, minusCor)
  
  ccleRNA <- ccleRNA[, names(selectCor)]
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    yOuterTrain <- thisFc[oi]
    
    innerIndex <- caret::createFolds(ytrain, k=10, list=TRUE, returnTrain = TRUE)
    
    innerMatrix <- sapply(innerIndex, function(ii){
      yInnerTrain <- ytrain[ii]
      yInnerTest <- ytrain[-ii]
      xInnerTrain <- xtrain[ii, ]
      xInnerTest <- xtrain[-ii, ]
      
      # glmnet package
      library(doParallel)
      cl <- makeCluster(10)
      registerDoParallel(cl)
      aGrid <- seq(0.1, 0.9, 0.05)
      search <- foreach(i=aGrid, .combine=rbind) %dopar% {
        cv <- glmnet::cv.glmnet(xInnerTrain, yInnerTrain, family="gaussian", nfold=10, alpha=i, type.measure = "deviance", parallel = TRUE)
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
      }
      cv3 <- search[search$cvm == min(search$cvm), ]
      stopCluster(cl)
      
      elastic_fit <- glmnet(xInnerTrain, yInnerTrain, family="gaussian", lambda = cv3$lambda.1se, alpha = cv3$alpha)
      lasso.coef <- predict(elastic_fit, type="coefficients", s=cv$lambda.lse)[,1]
      # lasso.coef <- lasso.coef[lasso.coef != 0]
      lasso.coef <- lasso.coef[2:length(lasso.coef)]
      
    })
    # return(innerMatrix)
    
  })
  
  return(outerList)
  
}



















