library(tidyverse)
library(purrr)

load("~/Data/mayo/project/project/rna2essential/DepQ4_select_gene_sgRNA_for_ccleLinnorm_prediction.RData")
m <- dplyr::filter(geneSummary2, maxRNA < 1, medianRNA == 0)
geneLeft <- dplyr::filter(geneSummary2, ! gene %in% m$gene)
ccleLinnorm <- ccleLinnorm[, geneLeft$gene]
ccleRNA <- ccleRNA[, geneLeft$gene]

setwd("/z/zhenyang/Data/mayo/project/project/rna2essential/BART_variable_selection")
bartList <- list.files(pattern="*_select.rds")

sapply(1:length(bartList), function(idx) {
  bartGene <- gsub("_BART_select.rds", "", bartList[1])
  bartFc <- theGeneFc[bartGene, ]
  bartFc <- (bartFc - min(bartFc))/(max(bartFc) - min(bartFc))   # normalize bartFc to range (0, 1)
  bartFeature <- readRDS(bartList[1])
  
  bartMatrix <- lapply(bartFeature, function(i) {
    lapply(i, function(j) j[, c("gene", "varCount")]) %>%
      purrr::reduce(full_join, by="gene")
  }) %>%
    purrr::reduce(full_join, by="gene") %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  bartCount <- apply(bartMatrix, 1, function(x) median(x, na.rm=TRUE)) %>% sort(decreasing = TRUE) %>% names()
  bartResult <- lapply(seq(8, 15, 1), function(y) {
    bartRNA <- ccleRNA[, bartCount[1:y]]
    bartOutput <- bart_predict2(bartRNA, bartFc, 20)
  })
  saveRDS(bartResult, paste0(bartGene, "_BART_tune.rds"))

})









# BART tune and predict

bart_predict2 <- function(ccleRNA, thisFc, folds=5, useCores=15) {
  
  
  outerIndex <- caret::createFolds(thisFc, k=folds, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    xOuterTest <- ccleRNA[-oi, ]
    yOuterTrain <- thisFc[oi]
    yOuterTest <- thisFc[-oi]
    
    combMatrix <- expand.grid(list(a=c(0.5, 0.7), k=c(2,3,4,5), sigq= c(0.9, 0.95, 0.99)))  %>% as.matrix()
    
    thisLoop <- lapply(1:nrow(combMatrix), function(x){
      a <- combMatrix[x, ][1]
      k <- combMatrix[x, ][2]
      sigquant <- combMatrix[x, ][3]
      
      dirichlet_select <- BART::mc.wbart(xOuterTrain, yOuterTrain, xOuterTest, a=a, cont=TRUE, sigdf=3, k=k, sigquant=sigquant, 
                                    ntree=200, ndpost=8000, nskip=2000, printevery=3000, mc.cores=useCores, nice=0)
      thisRMSE <- sqrt(mean((dirichlet_select$yhat.test.mean - yOuterTest)^2))
      thisCor <- cor(dirichlet_select$yhat.test.mean, yOuterTest, method="spearman")
      # rm(dirichlet_select)
      result <- list(thisRMSE, thisCor, dirichlet_select)
      })
    
    return(thisLoop)
  })
  return(outerList)
  
}


bart_predict <- function(ccleRNA, thisFc, folds=5, useCores=20) {
  
  outerIndex <- caret::createFolds(thisFc, k=folds, list=TRUE, returnTrain = TRUE)
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    xOuterTest <- ccleRNA[-oi, ]
    yOuterTrain <- thisFc[oi]
    yOuterTest <- thisFc[-oi]
    

    combMatrix <- expand.grid(list(a=c(0.5, 0.7), k=c(2,3,4,5), sigq= c(0.9, 0.95, 0.99)))  %>% as.matrix()
    thisLoop <- parallel::mclapply(1:nrow(combMatrix), function(x){
      a <- combMatrix[x, ][1]
      k <- combMatrix[x, ][2]
      sigquant <- combMatrix[x, ][3]
      dirichlet_select <- BART::wbart(xOuterTrain, yOuterTrain, xOuterTest, a=a, cont=TRUE, sigdf=3, k=k, sigquant=sigquant, 
                                      ntree=200, ndpost=5000, nskip=3000, printevery=1000)
      output <- list(dirichlet_select$varprob, dirichlet_select$varcount, dirichlet_select$varcount.mean, dirichlet_select$varprob.mean)
      
      thisRMSE <- sqrt(mean((dirichlet_select$yhat.test.mean - yOuterTest)^2))
      thisCor <- cor(dirichlet_select$yhat.test.mean, yOuterTest, method="spearman")
      rm(dirichlet_select)
      result <- list(rmse=thisRMSE, spearman=thisCor, modelOutput=output, parameter=c(a=a,k=k,sigquant=sigquant))
    }, mc.cores=useCores)
    # loopCor <- sapply(thisLoop, function(x) x[[2]])
    # theIdx <- which(loopCor == max(loopCor))
    #theSelect <- thisLoop[[theIdx]]
    
  })
  return(outerList)
  
}





### test computing time

dirichlet_select <- BART::mc.wbart(xOuterTrain, yOuterTrain, xOuterTest, cont=TRUE, sigdf=3, sigquant=0.9, k=3, 
                                                 ntree=200, ndpost=10000, nskip=2000, printevery=1000, mc.cores = 30, nice = 0)

dirichlet_select2 <- BART::wbart(xOuterTrain, yOuterTrain, xOuterTest,cont=TRUE, sigdf=3, sigquant=0.9, k=3, 
                                   ntree=200, ndpost=8000, nskip=2000, printevery=1000)




# caret-xgboost tune and predict
caret_xgb <- function(ccleRNA, thisFc, nthreads=28) {
  
  outerIndex <- caret::createFolds(thisFc, k=10, list=TRUE, returnTrain = TRUE)
  
  outerList <- lapply(outerIndex, function(oi){
    
    xOuterTrain <- ccleRNA[oi, ]
    yOuterTrain <- thisFc[oi]
    
    cv.ctrl <- caret::trainControl(method = "cv", repeats = 1, number = 10, 
                            allowParallel=TRUE)
    
    xgb.grid <- expand.grid(nrounds = 1000,
                            eta = c(0.01,0.05,0.1, 0.2),
                            max_depth = c(6,8,10,14,16,18,20),
                            gamma=0.2, 
                            colsample_bytree=1,
                            min_child_weight=1,
                            subsample=1)
    set.seed(45)
    xgb_tune <- caret::train(x=xOuterTrain,
                      y=yOuterTrain,
                      method="xgbTree",
                      trControl=cv.ctrl,
                      tuneGrid=xgb.grid,
                      verbose=T,
                      metric="RMSE",
                      nthread =nthreads)
    
  })
  
  return(outerList)
  
}

# xgboost tune and predict
xgb_predict <- function(xgbRNA, xgbFc, nthreads=28) {
  
  
  outerIndex <- caret::createFolds(xgbFc, k=10, list=TRUE, returnTrain = TRUE) 
  thisIdx <- outerIndex[[1]]
  xTrain <- xgbRNA[thisIdx, ]
  xTest <- xgbRNA[- thisIdx, ]
  yTrain <- xgbFc[thisIdx]
  yTest <- xgbFc[- thisIdx]
  
  eta <- 0.01
  colsample_bytree <- c(0.5, 0.7, 0.8, 0.9, 1)
  max_depth <- c(4, 6, 8, 10, 12)
  subsample=c(0.25, 0.5, 0.75, 1)
  gamma <- c(0, 0.1, 0.2, 0.3)
  
  
  dtrain <- xgboost::xgb.DMatrix(data=xTrain, label=as.numeric(yTrain))
  dtest <- xgboost::xgb.DMatrix(data=xTest, label=as.numeric(yTest))
  
  thisGrid <- expand.grid(max_depth, gamma, subsample, colsample_bytree) %>% setNames(c("max_depth", "gamma", "subsample", "colsample_bytree"))
 
  for(i in 1:nrow(thisGrid)) {
    params <- list(max_depth=thisGrid[i, 1], gamma=thisGrid[i, 2], subsample=thisGrid[i, 3], colsample_bytree=thisGrid[i, 4], nthread=20)
    xgb <- xgboost(params=params, data=dtrain, nrounds=3000, 
                      nfold=5, prediction=TRUE, metrics=list("rmse"), 
                      stratified=TRUE, print_every_n=500, 
                      early_stopping_rounds=50, objective="reg:linear", eta=0.01)
    conv_gamma[,i] = xgb$evaluation_log$train_rmse
    pred_gamma[,i] = predict(xgb, xTest)
  }
  
 
  
  watchlist <- list(train=dtrain, test=dtest)
  
      
      
      
      
      xgbSelect <- xgboost::xgb.train(data=dtrain, watchlist=watchlist, eta=eta, max.depth=depth, nrounds=5000, early_stopping_rounds = 50, print_every_n = 500, nthread=nthreads, verbose=TRUE)
      varGain <- xgboost::xgb.importance(model=xgbSelect) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("Feature") %>%
        as.matrix() %>%
        .[, 1]
      

}



# rrf
rrf_mlr_predict <- function(rrfRNA, rrfFc, useCPUs=10) {
  
  outerIndex <- caret::createFolds(rrfFc, k=10, list=TRUE, returnTrain = TRUE) 
  thisIdx <- outerIndex[[1]]
  xTrain <- rrfRNA[thisIdx, ]
  xTest <- rrfRNA[- thisIdx, ]
  yTrain <- rrfFc[thisIdx]
  yTest <- rrfFc[- thisIdx]
  
  
  trainDF <- data.frame(xTrain, response = yTrain)
  train_task <- makeRegrTask(data = trainDF, target="response")
  testDF <- data.frame(xTest, response = yTest)
  test_task <- makeRegrTask(data = testDF, target = 'response')
  
  cv_folds <- makeResampleDesc("CV", iters = 10)
  model <- makeLearner("regr.RRF", predict.type = "response")
  
  model_params <- makeParamSet(
    makeIntegerParam(id = "mtry", lower = 1, upper = 5), 
    makeDiscreteParam(id = "ntree", values = c(50, 100, 200, 300)),
    makeLogicalLearnerParam(id = "importance", default = TRUE),
    makeIntegerLearnerParam(id = "nodesize", lower = 5L, upper=15),
    makeIntegerLearnerParam(id = "maxnodes", lower = 5L, upper=15)
  )
  
  parallelMap::parallelStartMulticore(useCPUs)
  tuned_model <- tuneParams(learner = model, 
                            task = train_task,
                            resampling = cv_folds,
                            measure = list(rmse, rsq, spearmanrho),
                            par.set = model_params,
                            control =  makeTuneControlGrid(),
                            show.info = FALSE)
  
  
  lrn <- setHyperPars(makeLearner("regr.RRF"), par.vals = tuned_model$x)
  thisModel <- train(lrn, train_task)
  thisPredict <- predict(thisModel, task = test_task)
  
  # calculate performance
  thisPerformance <- performance(thisPredict, measures = list(rmse, rsq, spearmanrho))
  
  # stop parallel
  parallelMap::parallelStop()
  
  return(list(tuneModel=tuned_model, predict=thisPredict, performance=thisPerformance, test=testDF))
  
}


# ranger tune and predict
ranger_mlr_train <- function(rangerRNA, rangerFc, useCPUs=10) {
  require(mlr)
  # outerIndex <- caret::createFolds(rangerFc, k=10, list=TRUE, returnTrain = TRUE) 
  # thisIdx <- outerIndex[[1]]
  # xTrain <- rangerRNA[thisIdx, ]
  # xTest <- rangerRNA[- thisIdx, ]
  # yTrain <- rangerFc[thisIdx]
  # yTest <- rangerFc[- thisIdx]
  
  #create tasks
  trainDF <- data.frame(rangerRNA, response=rangerFc)
  train_task <- mlr::makeRegrTask(data = trainDF, target = 'response')
  
  # testDF <- data.frame(xTest, response=yTest)
  # test_task <- makeRegrTask(data = testDF, target = 'response')
  
  #create learner
  lrn = mlr::makeLearner('regr.ranger')
  
  cv_folds <- mlr::makeResampleDesc("CV", iters=10)
  model <- mlr::makeLearner("regr.ranger", predict.type="se", keep.inbag=TRUE)
  
  model_params <- makeParamSet(
    makeDiscreteParam("num.trees", values=2000), 
    makeIntegerParam("mtry", lower=1, upper=6),
    makeIntegerParam("max.depth", lower=6, upper=15), 
    makeDiscreteLearnerParam(id = "importance", values = "impurity"),
    makeLogicalLearnerParam(id = "replace", default = TRUE),
    makeLogicalLearnerParam(id = "write.forest", default= TRUE)
    #makeLogicalLearnerParam(id = "keep.inbag", default = TRUE, tunable = FALSE)
    #makeIntegerParam("max.depth", lower=6, upper=15)
  )
  
  parallelMap::parallelStartMulticore(useCPUs)
  tuned_model <- tuneParams(learner = model, 
                            task = train_task,
                            resampling = cv_folds,
                            measure = list(rmse, rsq, spearmanrho),
                            par.set = model_params,
                            control =  makeTuneControlGrid(),
                            show.info = FALSE)
  
  
  lrn <- setHyperPars(makeLearner('regr.ranger'), par.vals=tuned_model$x)
  thisModel <- train(lrn, train_task)
  
  trainPredict <- mlr::predict(thisModel, tast=train_task)
  trainPerformance <- mlr::performance(trainPredict, measures = list(rmse, rsq, spearmanrho))
  
  parallelMap::parallelStop()
  return(list(trainModel=thisModel, trainPredict=trainPredict, trainPerformance=trainPerformance, trainDF=trainDF))
  # return(thisModel)

  
  # thisPredict <- predict(thisModel, task=test_task)
  # thisPerformance <- performance(thisPredict, measures = list(rmse, rsq, spearmanrho))
  
  # stop parallel
  # parallelMap::parallelStop()
  
  # return(list(tuneModel=tuned_model, predict=thisPredict, performance=thisPerformance, test=testDF))
}


ranger_mlr_predict <- function(testRNA, testFc) {
  require(mlr)
  testDF <- data.frame(testRNA, response=testFc)
  test_task <- makeRegrTask(data = testDF, target = 'response')
  
  thisPredict <- predict(thisModel, task=test_task)
  
  # calculate performance
  thisPerformance <- performance(thisPredict, measures = list(rmse, rsq, spearmanrho))
  
  # stop parallel
  parallelMap::parallelStop()
  
  return(list(testPredict=thisPredict, testPperformance=thisPerformance, testDF=testDF))
}



# xgboost tune and predict
xgb_mlr_train <- function(xgbRNA, xgbFc, useCPUs=20) {
  
  #outerIndex <- caret::createFolds(xgbFc, k=10, list=TRUE, returnTrain = TRUE) 
  #thisIdx <- outerIndex[[1]]
  #xTrain <- xgbRNA[thisIdx, ]
  #xTest <- xgbRNA[- thisIdx, ]
  #yTrain <- xgbFc[thisIdx]
  #yTest <- xgbFc[- thisIdx]
  
  #create tasks
  trainDF <- data.frame(xgbRNA, response=xgbFc)
  train_task <- makeRegrTask(data = trainDF, target = 'response')
  # testDF <- data.frame(xTest, response=yTest)
  # test_task <- makeRegrTask(data = testDF, target = 'response')
  
  #create learner
  lrn = makeLearner('regr.xgboost', nthread=useCPUs, print_every_n=2000, nrounds=5000)
  
  #set parameter space
  params <- makeParamSet(
    makeDiscreteParam(id = "max_depth", values = c(4, 6, 8, 10, 12)),
    makeDiscreteParam(id = "min_child_weight", values = 1), # c(1, 3, 5, 7)) 
    makeDiscreteParam(id = "gamma",  values = c(0, 0.1, 0.2, 0.3)),
    makeDiscreteParam(id = "eta", values = 0.05),
    makeDiscreteParam(id = "colsample_bytree", values = c(0.5, 0.7, 0.8, 0.9, 1))
    #makeDiscreteParam(id = "print_every_n", values = 2000, tunable = FALSE)
  )
  
  #set resampling strategy
  rdesc <- makeResampleDesc("CV", iters = 10L)
  
  #search strategy
  ctrl <- makeTuneControlRandom(maxit = 100)
  
  #set parallel backend
  #if(Sys.info()['sysname'] == "Linux") {
  #  parallelMap::parallelStartMulticore(cpus = useCPUs, show.info = FALSE)
  #} else parallelMap::parallelStartSocket(cpus = useCPUs, show.info = FALSE)
  
  tuned_model <- tuneParams(learner = lrn,
                     task = train_task, 
                     resampling = rdesc, 
                     measures = list(rmse, rsq, spearmanrho), 
                     par.set = params, 
                     control =  makeTuneControlRandom(maxit=5L), 
                     show.info = FALSE)
  
  
  
  lrn <- setHyperPars(makeLearner('regr.xgboost'), par.vals=tuned_model$x)
  thisModel <- train(lrn, train_task)
  
  trainPredict <- mlr::predict(thisModel, tast=train_task)
  trainPerformance <- mlr::performance(trainPredict, measures = list(rmse, rsq, spearmanrho))
  
  return(list(trainModel=thisModel, trainPredict=trainPredict, trainPerformance=trainPerformance, trainDF=trainDF))

}


xgb_mlr_predict <- function(xgb_mlr_model, testRNA, testFc) {
  testDF <- data.frame(testRNA, response=testFc)
  test_task <- mlr::makeRegrTask(data = testDF, target = 'response')
  
  thisPredict <- mlr::predict(xgb_mlr_model, task=test_task)
  
  thisPerformance <- mlr::performance(thisPredict, measures = list(rmse, rsq, spearmanrho))
  
  
  return(list(testPredict=thisPredict, testPperformance=thisPerformance, testDF=testDF))
  
}



### bartMachine tune and predict
bartMachine_train <- function(bartMRNA, bartMFc, folds=5, useCPUs=20) {
  # create train and test data
  # outerIndex <- caret::createFolds(bartFc, k=10, list=TRUE, returnTrain = TRUE) 
  # thisIdx <- outerIndex[[1]]
  # xTrain <- bartRNA[thisIdx, ] %>% as.data.frame()
  # xTest <- bartRNA[- thisIdx, ] %>% as.data.frame()
  # yTrain <- bartFc[thisIdx] 
  # yTest <- bartFc[- thisIdx]
  
  
  # set parameters
  nu_q_cvs <- list(c(2,0.8), c(2, 0.9), c(2, 0.99), c(3, 0.8), c(3, 0.9), c(3, 0.99), c(4, 0.8), c(4, 0.9), c(4, 0.99), c(5, 0.8), c(5, 0.9), c(5, 0.99))
  
  bartMachine::set_bart_machine_num_cores(useCPUs)
  thisBartM_cv <- bartMachine::bartMachineCV(data.frame(bartMRNA), bartMFc, num_tree_cvs=c(150, 200), k_cvs=c(2,3,4),
                                             nu_q_cvs=nu_q_cvs, k_folds=folds, verbose=FALSE)
  
  trainPredict <- bartMachine::bart_predict_for_test_data(thisBartM_cv, Xtest=data.frame(bartMRNA), ytest=bartMFc)
  # thisCor <- cor(thisPredict$y_hat, yTest, method="spearman")
  
  return(list(trainModel=thisBartM_cv, trainPredict=trainPredict, trainDF=data.frame(bartMRNA, response=bartMFc)))
  # return(list(cvModel=thisBartM_cv, predict=thisPredict, spearmanRho=thisCor, test=data.frame(xTest, response=yTest)))
  
}

### bartMachine tune and predict
bartMachine_predict <- function(bartMmodel, testRNA, testFc) {
  
  testPredict <- bartMachine::bart_predict_for_test_data(bartMmodel, Xtest=data.frame(testRNA), ytest=testFc)
  # plot_y_vs_yhat(bartMmodel, Xtest=data.frame(testRNA), ytest=testFc)
  
  return(list(testPredict=testPredict, testDF=data.frame(testRNA, response=testFc)))
  
}



# bartMachine_mlr train 
bartMachine_mlr_train <- function(bartRNA, bartFc, useCPUs=10) {
  require(mlr)
  # create train and test data
  # outerIndex <- caret::createFolds(bartFc, k=10, list=TRUE, returnTrain = TRUE) 
  # thisIdx <- outerIndex[[1]]
  # xTrain <- bartRNA[thisIdx, ]
  # xTest <- bartRNA[- thisIdx, ]
  # yTrain <- bartFc[thisIdx]
  # yTest <- bartFc[- thisIdx]
  
  
  trainDF <- data.frame(bartRNA, response=bartFc) 
  # testDF <- data.frame(xTest, response=yTest) 
  train_task <- makeRegrTask(data = trainDF, target="response")
  # test_task <- makeRegrTask(data = testDF, target="response")
  cv_folds <- makeResampleDesc("CV", iters=10)
  model <- makeLearner("regr.bartMachine", predict.type="response")
  
  model_params <- makeParamSet(
    makeDiscreteParam(id = "q", values = c(0.8, 0.9, 0.99)), 
    makeDiscreteParam(id = "num_trees", values = c(100, 150, 200)),
    makeDiscreteParam(id = "k", values = c(2, 3, 4, 5)),
    makeDiscreteParam(id = "nu", values = c(2, 3, 4, 5))
  )
  
  # parallel processing
  parallelMap::parallelStartMulticore(useCPUs)
  
  # tune parameters
  tuned_model <- tuneParams(learner = model, 
                            task = train_task,
                            resampling = cv_folds,
                            measure = list(rmse, spearmanrho, mse),
                            par.set = model_params,
                            control =  makeTuneControlGrid(),
                            show.info = FALSE)
  
  
  lrn <- setHyperPars(makeLearner("regr.bartMachine"), par.vals=tuned_model$x)
  thisModel <- train(lrn, regr.task)
  
  trainPredict <- mlr::predict(thisModel, tast=train_task)
  trainPerformance <- mlr::performance(trainPredict, measures = list(rmse, rsq, spearmanrho))
  # stop parallel
  parallelMap::parallelStop()
  
  return(list(trainModel=thisModel, trainPredict=trainPredict, trainPerformance=trainPerformance, trainDF=trainDF))
  
  
  # return(tuned_model)

}

# bartMachine_mlr predict
bartMachine_mlr_predcit <- function(bartMachineModel, testRNA, testFc) {
  require(mlr)
  testDF <- data.frame(testRNA, response=testFc) 
  test_task <- mlr::makeRegrTask(data = testDF, target="response")
  thisPredict <- mlr::predict(bartMachineModel, task=test_task)
  
  # calculate performance
  thisPerformance <- mlr::performance(thisPredict, measures = list(rmse, rsq, spearmanrho))

  return(list(testPredict=thisPredict, testPperformance=thisPerformance, testDF=testDF))
}

# random forest train

# random forest predict