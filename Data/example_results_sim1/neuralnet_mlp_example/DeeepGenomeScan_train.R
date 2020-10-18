"train" <-
  function(x, ...){
    UseMethod("train")
  }

#' @rdname train
#' @importFrom stats predict
#' @importFrom utils object.size flush.console
#' @importFrom withr with_seed
#' @export
train.default <- function(x, y,
                          method = "rf",
                          preProcess = NULL,
                          ...,
                          weights = NULL,
                          metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                          maximize = ifelse(metric %in% c("RMSE", "logLoss", "MAE"), FALSE, TRUE),
                          trControl = trainControl(),
                          tuneGrid = NULL,
                          tuneLength = ifelse(trControl$method == "none", 1, 3)) {
  startTime <- proc.time()
  
  ## get a seed before packages are loaded or recipes are processed
  rs_seed <- sample.int(.Machine$integer.max, 1L)
  
  if(is.null(colnames(x)))
    stop("Please use column names for `x`", call. = FALSE)
  
  if(is.character(y)) y <- as.factor(y)
  
  if( !is.numeric(y) & !is.factor(y) ){
    msg <- paste("Please make sure that the outcome column is a factor or numeric.",
                 "The class(es) of the column:",
                 paste0("'", class(y), "'", collapse = ", "))
    
    stop(msg, call. = FALSE )
  }
  
  if(is.list(method)) {
    minNames <- c("library", "type", "parameters", "grid",
                  "fit", "predict", "prob")
    nameCheck <- minNames %in% names(method)
    if(!all(nameCheck)) stop(paste("some required components are missing:",
                                   paste(minNames[!nameCheck], collapse = ", ")),
                             call. = FALSE)
    models <- method
    method <- "custom"
  } else {
    models <- getModelInfo(method, regex = FALSE)[[1]]
    if (length(models) == 0)
      stop(paste("Model", method, "is not in caret's built-in library"), call. = FALSE)
  }
  checkInstall(models$library)
  for(i in seq(along = models$library)) do.call("requireNamespaceQuietStop", list(package = models$library[i]))
  if(any(names(models) == "check") && is.function(models$check)) {
    software_check <- models$check(models$library)
  }
  
  
  paramNames <- as.character(models$parameters$parameter)
  
  funcCall <- match.call(expand.dots = TRUE)
  modelType <- get_model_type(y)
  if(!(modelType %in% models$type)) stop(paste("wrong model type for", tolower(modelType)), call. = FALSE)
  
  if(grepl("^svm", method) & grepl("String$", method)) {
    if(is.vector(x) && is.character(x)) {
      stop("'x' should be a character matrix with a single column for string kernel methods", call. = FALSE)
    }
    if(is.matrix(x) && is.numeric(x)) {
      stop("'x' should be a character matrix with a single column for string kernel methods", call. = FALSE)
    }
    if(is.data.frame(x)) {
      stop("'x' should be a character matrix with a single column for string kernel methods", call. = FALSE)
    }
  }
  
  if(modelType == "Regression" & length(unique(y)) == 2)
    warning(paste("You are trying to do regression and your outcome only has",
                  "two possible values Are you trying to do classification?",
                  "If so, use a 2 level factor as your outcome column."))
  
  if(modelType != "Classification" & !is.null(trControl$sampling))
    stop("sampling methods are only implemented for classification problems", call. = FALSE)
  if(!is.null(trControl$sampling)) {
    trControl$sampling <- parse_sampling(trControl$sampling)
  }
  
  if(any(class(x) == "data.table")) x <- as.data.frame(x, stringsAsFactors = TRUE)
  check_dims(x = x, y = y)
  n <- if(class(y)[1] == "Surv") nrow(y) else length(y)
  
  ## TODO add check method and execute here
  
  ## Some models that use RWeka start multiple threads and this conflicts with multicore:
  parallel_check("RWeka", models)
  parallel_check("keras", models)
  
  if(!is.null(preProcess) && !(all(names(preProcess) %in% ppMethods)))
    stop(paste('pre-processing methods are limited to:', paste(ppMethods, collapse = ", ")), call. = FALSE)
  if(modelType == "Classification") {
    ## We should get and save the class labels to ensure that predictions are coerced
    ## to factors that have the same levels as the original data. This is especially
    ## important with multiclass systems where one or more classes have low sample sizes
    ## relative to the others
    classLevels <- levels(y)
    attributes(classLevels) <- list(ordered = is.ordered(y))
    xtab <- table(y)
    if(any(xtab == 0)) {
      xtab_msg <- paste("'", names(xtab)[xtab == 0], "'", collapse = ", ", sep = "")
      stop(paste("One or more factor levels in the outcome has no data:", xtab_msg), call. = FALSE)
    }
    
    if(trControl$classProbs && any(classLevels != make.names(classLevels))) {
      stop(paste("At least one of the class levels is not a valid R variable name;",
                 "This will cause errors when class probabilities are generated because",
                 "the variables names will be converted to ",
                 paste(make.names(classLevels), collapse = ", "),
                 ". Please use factor levels that can be used as valid R variable names",
                 " (see ?make.names for help)."), call. = FALSE)
    }
    
    if(metric %in% c("RMSE", "Rsquared"))
      stop(paste("Metric", metric, "not applicable for classification models"), call. = FALSE)
    if(!trControl$classProbs && metric == "ROC")
      stop(paste("Class probabilities are needed to score models using the",
                 "area under the ROC curve. Set `classProbs = TRUE`",
                 "in the trainControl() function."), call. = FALSE)
    
    if(trControl$classProbs) {
      if(!is.function(models$prob)) {
        warning("Class probabilities were requested for a model that does not implement them")
        trControl$classProbs <- FALSE
      }
    }
  } else {
    if(metric %in% c("Accuracy", "Kappa"))
      stop(paste("Metric", metric, "not applicable for regression models"), call. = FALSE)
    classLevels <- NA
    if(trControl$classProbs) {
      warning("cannnot compute class probabilities for regression")
      trControl$classProbs <- FALSE
    }
  }
  
  
  if(trControl$method == "oob" & is.null(models$oob))
    stop("Out of bag estimates are not implemented for this model", call. = FALSE)
  
  ## If they don't exist, make the data partitions for the resampling iterations.
  trControl <- withr::with_seed(
    rs_seed,
    make_resamples(trControl, outcome = y)
  )
  
  if(is.logical(trControl$savePredictions)) {
    trControl$savePredictions <- if(trControl$savePredictions) "all" else "none"
  } else {
    if(!(trControl$savePredictions %in% c("all", "final", "none")))
      stop('`savePredictions` should be either logical or "all", "final" or "none"', call. = FALSE)
  }
  
  ## Gather all the pre-processing info. We will need it to pass into the grid creation
  ## code so that there is a concordance between the data used for modeling and grid creation
  if(!is.null(preProcess)) {
    ppOpt <- list(options = preProcess)
    if(length(trControl$preProcOptions) > 0) ppOpt <- c(ppOpt,trControl$preProcOptions)
  } else ppOpt <- NULL
  
  ## If no default training grid is specified, get one. We have to pass in the formula
  ## and data for some models (rpart, pam, etc - see manual for more details)
  if(is.null(tuneGrid)) {
    if(!is.null(ppOpt) && length(models$parameters$parameter) > 1 &&
       as.character(models$parameters$parameter) != "parameter") {
      pp <- list(method = ppOpt$options)
      if("ica" %in% pp$method) pp$n.comp <- ppOpt$ICAcomp
      if("pca" %in% pp$method) pp$thresh <- ppOpt$thresh
      if("knnImpute" %in% pp$method) pp$k <- ppOpt$k
      pp$x <- x
      ppObj <- do.call("preProcess", pp)
      tuneGrid <- models$grid(x = predict(ppObj, x),
                              y = y,
                              len = tuneLength,
                              search = trControl$search)
      rm(ppObj, pp)
    } else {
      tuneGrid <- models$grid(x = x, y = y, len = tuneLength, search = trControl$search)
      if (trControl$search != "grid" && tuneLength < nrow(tuneGrid))
        tuneGrid <- tuneGrid[1:tuneLength,,drop = FALSE]
    }
  }
  
  ## Remove duplicates from grid that can occur with random sampling and discrete model parameters
  tuneGrid <- tuneGrid[!duplicated(tuneGrid), , drop = FALSE]
  
  ## Check to make sure that there are tuning parameters in some cases
  if(grepl("adaptive", trControl$method) & nrow(tuneGrid) == 1) {
    stop(paste("For adaptive resampling, there needs to be more than one",
               "tuning parameter for evaluation"), call. = FALSE)
  }
  
  dotNames <- hasDots(tuneGrid, models)
  if(dotNames) colnames(tuneGrid) <- gsub("^\\.", "", colnames(tuneGrid))
  ## Check tuning parameter names
  tuneNames <- as.character(models$parameters$parameter)
  if (!all(colnames(tuneGrid) %in% tuneNames)) {
    badNames <- colnames(tuneGrid)[!(colnames(tuneGrid) %in% tuneNames)]
    stop(paste("The tuning parameter grid should have columns",
               paste(tuneNames, collapse = ", ", sep = "")), call. = FALSE)
  }
  
  goodNames <- all.equal(sort(tuneNames), sort(names(tuneGrid)))
  
  if(!is.logical(goodNames) || !goodNames) {
    stop(paste("The tuning parameter grid should have columns",
               paste(tuneNames, collapse = ", ", sep = "")), call. = FALSE)
  }
  
  if(trControl$method == "none" && nrow(tuneGrid) != 1)
    stop("Only one model should be specified in tuneGrid with no resampling", call. = FALSE)
  
  
  ## In case prediction bounds are used, compute the limits. For now,
  ## store these in the control object since that gets passed everywhere
  trControl$yLimits <- if(is.numeric(y)) get_range(y) else NULL
  
  
  if(trControl$method != "none") {
    ##------------------------------------------------------------------------------------------------------------------------------------------------------#
    
    ## For each tuning parameter combination, we will loop over them, fit models and generate predictions.
    ## We only save the predictions at this point, not the models (and in the case of method = "oob" we
    ## only save the prediction summaries at this stage.
    
    ## trainInfo will hold the information about how we should loop to train the model and what types
    ## of parameters are used.
    
    ## There are two types of methods to build the models: "basic" means that each tuning parameter
    ## combination requires it's own model fit and "seq" where a single model fit can be used to
    ## get predictions for multiple tuning parameters.
    
    ## The tuneScheme() function is in miscr.R and it helps define the following:
    ##   - A data frame called "loop" with columns for parameters and a row for each model to be fit.
    ##     For "basic" models, this is the same as the tuning grid. For "seq" models, it is only
    ##     the subset of parameters that need to be fit
    ##   - A list called "submodels". If "basic", it is NULL. For "seq" models, it is a list. Each list
    ##     item is a data frame of the parameters that need to be varied for the corresponding row of
    ##     the loop object.
    ##
    ## For example, for a gbm model, our tuning grid might be:
    ##    .interaction.depth .n.trees .shrinkage
    ##                     1       50        0.1
    ##                     1      100        0.1
    ##                     2       50        0.1
    ##                     2      100        0.1
    ##                     2      150        0.1
    ##
    ## For this example:
    ##
    ##   loop:
    ##   .interaction.depth .shrinkage .n.trees
    ##                    1        0.1      100
    ##                    2        0.1      150
    ##
    ##   submodels:
    ##   [[1]]
    ##     .n.trees
    ##           50
    ##
    ##   [[2]]
    ##     .n.trees
    ##           50
    ##          100
    ##
    ## A simplified version of predictionFunction() would have the following gbm section:
    ##
    ##     # First get the predictions with the value of n.trees as given in the current
    ##     # row of loop
    ##     out <- predict(modelFit,
    ##                    newdata,
    ##                    type = "response",
    ##                    n.trees = modelFit$tuneValue$.n.trees)
    ##
    ##     # param is the current value of submodels. In normal prediction mode (i.e
    ##     # when using predict.train), param = NULL. When called within train()
    ##     # with this model, it will have the other values for n.trees.
    ##     # In this case, the output of the function is a list of predictions
    ##     # These values are deconvoluted in workerTasks() in misc.R
    ##     if(!is.null(param))
    ##       {
    ##         tmp <- vector(mode = "list", length = nrow(param) + 1)
    ##         tmp[[1]] <- out
    ##
    ##         for(j in seq(along = param$.n.trees))
    ##           {
    ##             tmp[[j]]  <- predict(modelFit,
    ##                                  newdata,
    ##                                  type = "response",
    ##                                  n.trees = param$.n.trees[j])
    ##           }
    ##         out <- tmp
    ##
    
    # paramCols <- paste(".", as.character(models$parameters$parameter), sep = "")
    
    if(is.function(models$loop) && nrow(tuneGrid) > 1){
      trainInfo <- models$loop(tuneGrid)
      if(!all(c("loop", "submodels") %in% names(trainInfo)))
        stop("The 'loop' function should produce a list with elements 'loop' and 'submodels'", call. = FALSE)
      lengths <- unlist(lapply(trainInfo$submodels, nrow))
      if(all(lengths == 0)) trainInfo$submodels <- NULL
    } else trainInfo <- list(loop = tuneGrid)
    
    
    num_rs <- if(trControl$method != "oob") length(trControl$index) else 1L
    if(trControl$method %in% c("boot632", "optimism_boot", "boot_all")) num_rs <- num_rs + 1L
    ## Set or check the seeds when needed
    if(is.null(trControl$seeds) || all(is.na(trControl$seeds)))  {
      seeds <- sample.int(n = 1000000L, size = num_rs * nrow(trainInfo$loop) + 1L)
      seeds <- lapply(seq(from = 1L, to = length(seeds), by = nrow(trainInfo$loop)),
                      function(x) { seeds[x:(x+nrow(trainInfo$loop)-1L)] })
      seeds[[num_rs + 1L]] <- seeds[[num_rs + 1L]][1L]
      trControl$seeds <- seeds
    } else {
      if(!(length(trControl$seeds) == 1 && is.na(trControl$seeds))) {
        ## check versus number of tasks
        numSeeds <- unlist(lapply(trControl$seeds, length))
        badSeed <- (length(trControl$seeds) < num_rs + 1L) ||
          (any(numSeeds[-length(numSeeds)] < nrow(trainInfo$loop))) ||
          (numSeeds[length(numSeeds)] < 1L)
        if(badSeed) stop(paste("Bad seeds: the seed object should be a list of length",
                               num_rs + 1, "with",
                               num_rs, "integer vectors of size",
                               nrow(trainInfo$loop), "and the last list element having at least a",
                               "single integer"), call. = FALSE)
        if(any(is.na(unlist(trControl$seeds)))) stop("At least one seed is missing (NA)", call. = FALSE)
      }
    }
    
    
    ## SURV TODO: modify defaultSummary for Surv objects
    if(trControl$method == "oob") {
      ## delay this test until later
      perfNames <- metric
    } else {
      ## run some data thru the summary function and see what we get
      testSummary <- evalSummaryFunction(y, wts = weights, ctrl = trControl,
                                         lev = classLevels, metric = metric,
                                         method = method)
      perfNames <- names(testSummary)
    }
    
    if(!(metric %in% perfNames)){
      oldMetric <- metric
      metric <- perfNames[1]
      warning(paste("The metric \"",
                    oldMetric,
                    "\" was not in ",
                    "the result set. ",
                    metric,
                    " will be used instead.",
                    sep = ""))
    }
    
    if(trControl$method == "oob"){
      tmp <- oobTrainWorkflow(x = x, y = y, wts = weights,
                              info = trainInfo, method = models,
                              ppOpts = preProcess, ctrl = trControl, lev = classLevels, ...)
      performance <- tmp
      perfNames <- colnames(performance)
      perfNames <- perfNames[!(perfNames %in% as.character(models$parameters$parameter))]
      if(!(metric %in% perfNames)){
        oldMetric <- metric
        metric <- perfNames[1]
        warning(paste("The metric \"",
                      oldMetric,
                      "\" was not in ",
                      "the result set. ",
                      metric,
                      " will be used instead.",
                      sep = ""))
      }
    } else {
      if(trControl$method == "LOOCV"){
        tmp <- looTrainWorkflow(x = x, y = y, wts = weights,
                                info = trainInfo, method = models,
                                ppOpts = preProcess, ctrl = trControl, lev = classLevels, ...)
        performance <- tmp$performance
      } else {
        if(!grepl("adapt", trControl$method)){
          tmp <- nominalTrainWorkflow(x = x, y = y, wts = weights,
                                      info = trainInfo, method = models,
                                      ppOpts = preProcess, ctrl = trControl, lev = classLevels, ...)
          performance <- tmp$performance
          resampleResults <- tmp$resample
        } else {
          tmp <- adaptiveWorkflow(x = x, y = y, wts = weights,
                                  info = trainInfo, method = models,
                                  ppOpts = preProcess,
                                  ctrl = trControl,
                                  lev = classLevels,
                                  metric = metric,
                                  maximize = maximize,
                                  ...)
          performance <- tmp$performance
          resampleResults <- tmp$resample
        }
      }
    }
    
    ## Remove extra indices
    trControl$indexExtra <- NULL
    
    ## TODO we used to give resampled results for LOO
    if(!(trControl$method %in% c("LOOCV", "oob"))) {
      if(modelType == "Classification" && length(grep("^\\cell", colnames(resampleResults))) > 0) {
        resampledCM <- resampleResults[, !(names(resampleResults) %in% perfNames)]
        resampleResults <- resampleResults[, -grep("^\\cell", colnames(resampleResults))]
        #colnames(resampledCM) <- gsub("^\\.", "", colnames(resampledCM))
      } else resampledCM <- NULL
    } else resampledCM <- NULL
    
    
    if(trControl$verboseIter)  {
      cat("Aggregating results\n")
      flush.console()
    }
    
    perfCols <- names(performance)
    perfCols <- perfCols[!(perfCols %in% paramNames)]
    
    if(all(is.na(performance[, metric]))) {
      cat(paste("Something is wrong; all the", metric, "metric values are missing:\n"))
      print(summary(performance[, perfCols[!grepl("SD$", perfCols)], drop = FALSE]))
      stop("Stopping", call. = FALSE)
    }
    
    ## Sort the tuning parameters from least complex to most complex
    if(!is.null(models$sort)) performance <- models$sort(performance)
    
    if(any(is.na(performance[, metric])))
      warning("missing values found in aggregated results")
    
    
    if(trControl$verboseIter && nrow(performance) > 1) {
      cat("Selecting tuning parameters\n")
      flush.console()
    }
    
    ## select the optimal set
    selectClass <- class(trControl$selectionFunction)[1]
    
    ## Select the "optimal" tuning parameter.
    if(grepl("adapt", trControl$method)) {
      perf_check <- subset(performance, .B == max(performance$.B))
    } else perf_check <- performance
    
    ## Make adaptive only look at parameters with B = max(B)
    if(selectClass == "function") {
      bestIter <- trControl$selectionFunction(x = perf_check,
                                              metric = metric,
                                              maximize = maximize)
    }
    else {
      if(trControl$selectionFunction == "oneSE") {
        bestIter <- oneSE(perf_check,
                          metric,
                          length(trControl$index),
                          maximize)
      } else {
        bestIter <- do.call(trControl$selectionFunction,
                            list(x = perf_check,
                                 metric = metric,
                                 maximize = maximize))
      }
    }
    
    if(is.na(bestIter) || length(bestIter) != 1) stop("final tuning parameters could not be determined", call. = FALSE)
    
    if(grepl("adapt", trControl$method)) {
      best_perf <- perf_check[bestIter,as.character(models$parameters$parameter),drop = FALSE]
      performance$order <- 1:nrow(performance)
      bestIter <- merge(performance, best_perf)$order
      performance$order <- NULL
    }
    
    
    ## Based on the optimality criterion, select the tuning parameter(s)
    bestTune <- performance[bestIter, paramNames, drop = FALSE]
  } else {
    bestTune <- tuneGrid
    performance <- evalSummaryFunction(y, wts = weights,
                                       ctrl = trControl,
                                       lev = classLevels,
                                       metric = metric,
                                       method = method)
    perfNames <- names(performance)
    performance <- as.data.frame(t(performance), stringsAsFactors = TRUE)
    performance <- cbind(performance, tuneGrid)
    performance <- performance[-1,,drop = FALSE]
    tmp <- resampledCM <- NULL
  }
  ## Save some or all of the resampling summary metrics
  if(!(trControl$method %in% c("LOOCV", "oob", "none"))) {
    byResample <- switch(trControl$returnResamp,
                         none = NULL,
                         all = {
                           out <- resampleResults
                           colnames(out) <- gsub("^\\.", "", colnames(out))
                           out
                         },
                         final = {
                           out <- merge(bestTune, resampleResults)
                           out <- out[,!(names(out) %in% names(tuneGrid)), drop = FALSE]
                           out
                         })
  } else {
    byResample <- NULL
  }
  
  # names(bestTune) <- paste(".", names(bestTune), sep = "")
  
  ## Reorder rows of performance
  orderList <- list()
  for(i in seq(along = paramNames)) orderList[[i]] <- performance[,paramNames[i]]
  
  performance <- performance[do.call("order", orderList),]
  
  if(trControl$verboseIter) {
    bestText <- paste(paste(names(bestTune), "=",
                            format(bestTune, digits = 3)),
                      collapse = ", ")
    if(nrow(performance) == 1) bestText <- "final model"
    cat("Fitting", bestText, "on full training set\n")
    flush.console()
  }
  
  ## Make the final model based on the tuning results
  
  indexFinal <- if(is.null(trControl$indexFinal)) seq(along = y) else trControl$indexFinal
  
  if(!(length(trControl$seeds) == 1 && is.na(trControl$seeds))) set.seed(trControl$seeds[[length(trControl$seeds)]][1])
  finalTime <- system.time(
    finalModel <- createModel(x = subset_x(x, indexFinal),
                              y = y[indexFinal],
                              wts = weights[indexFinal],
                              method = models,
                              tuneValue = bestTune,
                              obsLevels = classLevels,
                              pp = ppOpt,
                              last = TRUE,
                              classProbs = trControl$classProbs,
                              sampling = trControl$sampling,
                              ...))
  
  if(trControl$trim && !is.null(models$trim)) {
    if(trControl$verboseIter) old_size <- object.size(finalModel$fit)
    finalModel$fit <- models$trim(finalModel$fit)
    if(trControl$verboseIter) {
      new_size <- object.size(finalModel$fit)
      reduction <- format(old_size - new_size, units = "Mb")
      if(reduction == "0 Mb") reduction <- "< 0 Mb"
      p_reduction <- (unclass(old_size) - unclass(new_size))/unclass(old_size)*100
      p_reduction <- if(p_reduction < 1) "< 1%" else paste0(round(p_reduction, 0), "%")
      cat("Final model footprint reduced by", reduction, "or", p_reduction, "\n")
    }
  }
  
  ## get pp info
  pp <- finalModel$preProc
  finalModel <- finalModel$fit
  
  ## Remove this and check for other places it is reference
  ## replaced by tuneValue
  if(method == "pls") finalModel$bestIter <- bestTune
  
  ## To use predict.train and automatically use the optimal lambda,
  ## we need to save it
  if(method == "glmnet") finalModel$lambdaOpt <- bestTune$lambda
  
  if(trControl$returnData) {
    outData <- if(!is.data.frame(x)) try(as.data.frame(x, stringsAsFactors = TRUE), silent = TRUE) else x
    if(inherits(outData, "try-error")) {
      warning("The training data could not be converted to a data frame for saving")
      outData <- NULL
    } else   {
      outData$.outcome <- y
      if(!is.null(weights)) outData$.weights <- weights
    }
  } else outData <- NULL
  
  if(trControl$savePredictions == "final")
    tmp$predictions <- merge(bestTune, tmp$predictions)
  
  endTime <- proc.time()
  times <- list(everything = endTime - startTime,
                final = finalTime)
  
  out <- structure(list(method = method,
                        modelInfo = models,
                        modelType = modelType,
                        results = performance,
                        pred = tmp$predictions,
                        bestTune = bestTune,
                        call = funcCall,
                        dots = list(...),
                        metric = metric,
                        control = trControl,
                        finalModel = finalModel,
                        preProcess = pp,
                        trainingData = outData,
                        resample = byResample,
                        resampledCM = resampledCM,
                        perfNames = perfNames,
                        maximize = maximize,
                        yLimits = trControl$yLimits,
                        times = times,
                        levels = classLevels),
                   class = "train")
  trControl$yLimits <- NULL
  
  if(trControl$timingSamps > 0) {
    pData <- x[sample(1:nrow(x), trControl$timingSamps, replace = TRUE),,drop = FALSE]
    out$times$prediction <- system.time(predict(out, pData))
  } else  out$times$prediction <- rep(NA, 3)
  out
  
}

#' @rdname train
#' @importFrom stats .getXlevels complete.cases contrasts model.frame model.matrix model.response model.weights na.fail
#' @export
train.formula <- function (form, data, ..., weights, subset, na.action = na.fail, contrasts = NULL)  {
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval.parent(m$data)))  m$data <- as.data.frame(data, stringsAsFactors = TRUE)
  m$... <- m$contrasts <- NULL
  
  check_na_conflict(match.call(expand.dots = TRUE))
  
  ## Look for missing `na.action` in call. To make the default (`na.fail`)
  ## recognizable by `eval.parent(m)`, we need to add it to the call
  ## object `m`
  
  if(!("na.action" %in% names(m))) m$na.action <- quote(na.fail)
  
  # do we need the double colon here?
  m[[1]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  if(nrow(m) < 1) stop("Every row has at least one missing value were found", call. = FALSE)
  Terms <- attr(m, "terms")
  x <- model.matrix(Terms, m, contrasts)
  cons <- attr(x, "contrast")
  int_flag <- grepl("(Intercept)", colnames(x))
  if (any(int_flag)) x <- x[, !int_flag, drop = FALSE]
  w <- as.vector(model.weights(m))
  y <- model.response(m)
  
  res <- train(x, y, weights = w, ...)
  res$terms <- Terms
  res$coefnames <- colnames(x)
  res$call <- match.call()
  res$na.action <- attr(m, "na.action")
  res$contrasts <- cons
  res$xlevels <- .getXlevels(Terms, m)
  if(!is.null(res$trainingData)) {
    ## We re-save the original data from the formula interface
    ## since it has not been converted to dummy variables.
    res$trainingData <- data[,all.vars(Terms), drop = FALSE]
    isY <- names(res$trainingData) %in% as.character(form[[2]])
    if(any(isY)) colnames(res$trainingData)[isY] <- ".outcome"
  }
  class(res) <- c("train", "train.formula")
  res
}

#' @rdname train
#' @importFrom withr with_seed
#' @export
train.recipe <- function(x,
                         data,
                         method = "rf",
                         ...,
                         metric = ifelse(is.factor(y_dat), "Accuracy", "RMSE"),
                         maximize = ifelse(metric %in% c("RMSE", "logLoss", "MAE"), FALSE, TRUE),
                         trControl = trainControl(),
                         tuneGrid = NULL,
                         tuneLength = ifelse(trControl$method == "none", 1, 3)) {
  startTime <- proc.time()
  
  ## get a seed before packages are loaded or recipes are processed
  rs_seed <- sample.int(.Machine$integer.max, 1L)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  preproc_dots(...)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(is.list(method)) {
    minNames <- c("library", "type", "parameters", "grid",
                  "fit", "predict", "prob")
    nameCheck <- minNames %in% names(method)
    if(!all(nameCheck)) stop(paste("some required components are missing:",
                                   paste(minNames[!nameCheck], collapse = ", ")),
                             call. = FALSE)
    models <- method
    method <- "custom"
  } else {
    models <- getModelInfo(method, regex = FALSE)[[1]]
    if (length(models) == 0)
      stop(paste("Model", method, "is not in caret's built-in library"), call. = FALSE)
  }
  checkInstall(models$library)
  for(i in seq(along = models$library))
    do.call("requireNamespace", list(package = models$library[i]))
  if(any(names(models) == "check") && is.function(models$check)) {
    software_check <- models$check(models$library)
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # prep and bake recipe on entire training set
  if(trControl$verboseIter)  {
    cat("Preparing recipe\n")
    flush.console()
  }
  
  trained_rec <- prep(x, training = data,
                      fresh = TRUE,
                      retain = TRUE,
                      verbose = FALSE,
                      strings_as_factors = TRUE)
  x_dat <- juice(trained_rec, all_predictors())
  y_dat <- juice(trained_rec, all_outcomes())
  if(ncol(y_dat) > 1)
    stop("`train` doesn't support multivariate outcomes")
  y_dat <- getElement(y_dat, names(y_dat))
  is_weight <- summary(trained_rec)$role == "case weight"
  if(any(is_weight)) {
    if(sum(is_weight) > 1)
      stop("Ony one column can be used as a case weight.")
    weights <- juice(trained_rec, has_role("case weight"))
    weights <- getElement(weights, names(weights))
  } else weights <- NULL
  
  is_perf <- summary(trained_rec)$role == "performance var"
  if(any(is_perf)) {
    perf_data <- juice(trained_rec, has_role("performance var"))
  } else perf_data <- NULL
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  paramNames <- as.character(models$parameters$parameter)
  
  funcCall <- match.call(expand.dots = TRUE)
  modelType <- get_model_type(y_dat)
  if(!(modelType %in% models$type))
    stop(paste("wrong model type for", tolower(modelType)), call. = FALSE)
  
  ## RECIPE the rec might produce character `x_dat` so convert if these
  ## models are used? These need to be re-though since no matrix results
  if(grepl("^svm", method) & grepl("String$", method)) {
    if(is.vector(x_dat) && is.character(x_dat)) {
      stop("'x_dat' should be a character matrix with a single column for string kernel methods",
           call. = FALSE)
    }
    if(is.matrix(x_dat) && is.numeric(x_dat)) {
      stop("'x_dat' should be a character matrix with a single column for string kernel methods",
           call. = FALSE)
    }
    if(is.data.frame(x_dat)) {
      stop("'x_dat' should be a character matrix with a single column for string kernel methods",
           call. = FALSE)
    }
  }
  
  if(modelType == "Regression" & length(unique(y_dat)) == 2)
    warning(paste("You are trying to do regression and your outcome only has",
                  "two possible values Are you trying to do classification?",
                  "If so, use a 2 level factor as your outcome column."))
  
  if(modelType != "Classification" & !is.null(trControl$sampling))
    stop("sampling methods are only implemented for classification problems",
         call. = FALSE)
  if(!is.null(trControl$sampling)) {
    trControl$sampling <- parse_sampling(trControl$sampling)
  }
  
  check_dims(x = x_dat, y = y_dat)
  n <- if(class(y_dat)[1] == "Surv") nrow(y_dat) else length(y_dat)
  
  ## Some models that use RWeka start multiple threads and this conflicts with multicore:
  parallel_check("RWeka", models)
  parallel_check("keras", models)
  
  if(modelType == "Classification") {
    ## We should get and save the class labels to ensure that predictions are coerced
    ## to factors that have the same levels as the original data. This is especially
    ## important with multiclass systems where one or more classes have low sample sizes
    ## relative to the others
    classLevels <- levels(y_dat)
    attributes(classLevels) <- list(ordered = is.ordered(y_dat))
    xtab <- table(y_dat)
    if(any(xtab == 0)) {
      xtab_msg <- paste("'", names(xtab)[xtab == 0], "'", collapse = ", ", sep = "")
      stop(paste("One or more factor levels in the outcome has no data:", xtab_msg),
           call. = FALSE)
    }
    
    if(trControl$classProbs && any(classLevels != make.names(classLevels))) {
      stop(paste("At least one of the class levels is not a valid R variable name;",
                 "This will cause errors when class probabilities are generated because",
                 "the variables names will be converted to ",
                 paste(make.names(classLevels), collapse = ", "),
                 ". Please use factor levels that can be used as valid R variable names",
                 " (see ?make.names for help)."), call. = FALSE)
    }
    
    if(metric %in% c("RMSE", "Rsquared"))
      stop(paste("Metric", metric, "not applicable for classification models"),
           call. = FALSE)
    if(!trControl$classProbs && metric == "ROC")
      stop(paste("Class probabilities are needed to score models using the",
                 "area under the ROC curve. Set `classProbs = TRUE`",
                 "in the trainControl() function."), call. = FALSE)
    
    if(trControl$classProbs) {
      if(!is.function(models$prob)) {
        warning("Class probabilities were requested for a model that does not implement them")
        trControl$classProbs <- FALSE
      }
    }
  } else {
    if(metric %in% c("Accuracy", "Kappa"))
      stop(paste("Metric", metric, "not applicable for regression models"),
           call. = FALSE)
    classLevels <- NA
    if(trControl$classProbs) {
      warning("cannnot compute class probabilities for regression")
      trControl$classProbs <- FALSE
    }
  }
  
  if(trControl$method == "oob" & is.null(models$oob))
    stop("Out of bag estimates are not implemented for this model",
         call. = FALSE)
  
  ## If they don't exist, make the data partitions for the resampling iterations.
  # Get outcomes from the _original_ data since that is what should be given to
  # the recipe
  y_orig_val <- trained_rec$var_info$variable[trained_rec$var_info$role == "outcome"]
  y_orig_val <- y_orig_val
  trControl <- withr::with_seed(
    rs_seed,
    make_resamples(trControl, outcome = data[[y_orig_val]])
  )
  
  if(is.logical(trControl$savePredictions)) {
    trControl$savePredictions <- if(trControl$savePredictions) "all" else "none"
  } else {
    if(!(trControl$savePredictions %in% c("all", "final", "none")))
      stop('`savePredictions` should be either logical or "all", "final" or "none"', call. = FALSE)
  }
  
  if(is.null(tuneGrid)) {
    tuneGrid <- models$grid(x = x_dat, y = y_dat, len = tuneLength, search = trControl$search)
    if (trControl$search != "grid" && tuneLength < nrow(tuneGrid))
      tuneGrid <- tuneGrid[1:tuneLength,,drop = FALSE]
  }
  
  ## Check to make sure that there are tuning parameters in some cases
  if(grepl("adaptive", trControl$method) & nrow(tuneGrid) == 1) {
    stop(paste("For adaptive resampling, there needs to be more than one",
               "tuning parameter for evaluation"), call. = FALSE)
  }
  
  dotNames <- hasDots(tuneGrid, models)
  if(dotNames) colnames(tuneGrid) <- gsub("^\\.", "", colnames(tuneGrid))
  ## Check tuning parameter names
  tuneNames <- as.character(models$parameters$parameter)
  goodNames <- all.equal(sort(tuneNames), sort(names(tuneGrid)))
  
  if(!is.logical(goodNames) || !goodNames) {
    stop(paste("The tuning parameter grid should have columns",
               paste(tuneNames, collapse = ", ", sep = "")), call. = FALSE)
  }
  
  if(trControl$method == "none" && nrow(tuneGrid) != 1)
    stop("Only one model should be specified in tuneGrid with no resampling", call. = FALSE)
  
  ## In case prediction bounds are used, compute the limits. For now,
  ## store these in the control object since that gets passed everywhere
  trControl$yLimits <- if(is.numeric(y_dat)) get_range(y_dat) else NULL
  
  if(trControl$method != "none") {
    
    if(is.function(models$loop) && nrow(tuneGrid) > 1){
      trainInfo <- models$loop(tuneGrid)
      if(!all(c("loop", "submodels") %in% names(trainInfo)))
        stop("The 'loop' function should produce a list with elements 'loop' and 'submodels'", call. = FALSE)
      lengths <- unlist(lapply(trainInfo$submodels, nrow))
      if(all(lengths == 0)) trainInfo$submodels <- NULL
    } else trainInfo <- list(loop = tuneGrid)
    
    num_rs <- if(trControl$method != "oob") length(trControl$index) else 1L
    if(trControl$method %in% c("boot632", "optimism_boot", "boot_all")) num_rs <- num_rs + 1L
    ## Set or check the seeds when needed
    if(is.null(trControl$seeds) || all(is.na(trControl$seeds)))  {
      seeds <- sample.int(n = 1000000L, size = num_rs * nrow(trainInfo$loop) + 1L)
      seeds <- lapply(seq(from = 1L, to = length(seeds), by = nrow(trainInfo$loop)),
                      function(x) { seeds[x:(x+nrow(trainInfo$loop)-1L)] })
      seeds[[num_rs + 1L]] <- seeds[[num_rs + 1L]][1L]
      trControl$seeds <- seeds
    } else {
      if(!(length(trControl$seeds) == 1 && is.na(trControl$seeds))) {
        ## check versus number of tasks
        numSeeds <- unlist(lapply(trControl$seeds, length))
        badSeed <- (length(trControl$seeds) < num_rs + 1L) ||
          (any(numSeeds[-length(numSeeds)] < nrow(trainInfo$loop))) ||
          (numSeeds[length(numSeeds)] < 1L)
        if(badSeed) stop(paste("Bad seeds: the seed object should be a list of length",
                               num_rs + 1, "with",
                               num_rs, "integer vectors of size",
                               nrow(trainInfo$loop), "and the last list element having at least a",
                               "single integer"), call. = FALSE)
        if(any(is.na(unlist(trControl$seeds)))) stop("At least one seed is missing (NA)", call. = FALSE)
      }
    }
    if(trControl$method == "oob") {
      ## delay this test until later
      perfNames <- metric
    } else {
      ## run some data thru the summary function and see what we get
      testSummary <- evalSummaryFunction(y_dat,
                                         perf = perf_data,
                                         wts = weights, ctrl = trControl,
                                         lev = classLevels, metric = metric,
                                         method = method)
      perfNames <- names(testSummary)
    }
    
    if(!(metric %in% perfNames)){
      oldMetric <- metric
      metric <- perfNames[1]
      warning(paste("The metric \"",
                    oldMetric,
                    "\" was not in ",
                    "the result set. ",
                    metric,
                    " will be used instead.",
                    sep = ""))
    }
    if(trControl$method == "oob"){
      tmp <- oob_train_rec(rec = x, dat = data,
                           info = trainInfo, method = models,
                           ctrl = trControl, lev = classLevels, ...)
      performance <- tmp
      perfNames <- colnames(performance)
      perfNames <- perfNames[!(perfNames %in% as.character(models$parameters$parameter))]
      if(!(metric %in% perfNames)){
        oldMetric <- metric
        metric <- perfNames[1]
        warning(paste("The metric \"",
                      oldMetric,
                      "\" was not in ",
                      "the result set. ",
                      metric,
                      " will be used instead.",
                      sep = ""))
      }
    } else {
      if(trControl$method == "LOOCV"){
        tmp <- loo_train_rec(rec = x, dat = data,
                             info = trainInfo, method = models,
                             ctrl = trControl, lev = classLevels, ...)
        performance <- tmp$performance
      } else {
        if(!grepl("adapt", trControl$method)){
          tmp <- train_rec(rec = x, dat = data,
                           info = trainInfo, method = models,
                           ctrl = trControl, lev = classLevels, ...)
          performance <- tmp$performance
          resampleResults <- tmp$resample
        } else {
          tmp <- train_adapt_rec(rec = x, dat = data,
                                 info = trainInfo,
                                 method = models,
                                 ctrl = trControl,
                                 lev = classLevels,
                                 metric = metric,
                                 maximize = maximize,
                                 ...)
          performance <- tmp$performance
          resampleResults <- tmp$resample
        }
      }
    }
    
    ## Remove extra indices
    trControl$indexExtra <- NULL
    
    if(!(trControl$method %in% c("LOOCV", "oob"))) {
      if(modelType == "Classification" && length(grep("^\\cell", colnames(resampleResults))) > 0) {
        resampledCM <- resampleResults[, !(names(resampleResults) %in% perfNames)]
        resampleResults <- resampleResults[, -grep("^\\cell", colnames(resampleResults))]
        #colnames(resampledCM) <- gsub("^\\.", "", colnames(resampledCM))
      } else resampledCM <- NULL
    } else resampledCM <- NULL
    
    
    if(trControl$verboseIter)  {
      cat("Aggregating results\n")
      flush.console()
    }
    
    perfCols <- names(performance)
    perfCols <- perfCols[!(perfCols %in% paramNames)]
    
    if(all(is.na(performance[, metric]))) {
      cat(paste("Something is wrong; all the", metric, "metric values are missing:\n"))
      print(summary(performance[, perfCols[!grepl("SD$", perfCols)], drop = FALSE]))
      stop("Stopping", call. = FALSE)
    }
    
    ## Sort the tuning parameters from least complex to most complex
    if(!is.null(models$sort)) performance <- models$sort(performance)
    
    if(any(is.na(performance[, metric])))
      warning("missing values found in aggregated results")
    
    
    if(trControl$verboseIter && nrow(performance) > 1) {
      cat("Selecting tuning parameters\n")
      flush.console()
    }
    
    ## select the optimal set
    selectClass <- class(trControl$selectionFunction)[1]
    
    ## Select the "optimal" tuning parameter.
    if(grepl("adapt", trControl$method)) {
      perf_check <- subset(performance, Num_Resamples == max(performance$Num_Resamples))
    } else perf_check <- performance
    
    ## Make adaptive only look at parameters with B = max(B)
    if(selectClass == "function") {
      bestIter <- trControl$selectionFunction(x = perf_check,
                                              metric = metric,
                                              maximize = maximize)
    }
    else {
      if(trControl$selectionFunction == "oneSE") {
        bestIter <- oneSE(perf_check,
                          metric,
                          length(trControl$index),
                          maximize)
      } else {
        bestIter <- do.call(trControl$selectionFunction,
                            list(x = perf_check,
                                 metric = metric,
                                 maximize = maximize))
      }
    }
    
    if(is.na(bestIter) || length(bestIter) != 1)
      stop("final tuning parameters could not be determined", call. = FALSE)
    
    if(grepl("adapt", trControl$method)) {
      best_perf <- perf_check[bestIter,as.character(models$parameters$parameter),drop = FALSE]
      performance$order <- 1:nrow(performance)
      bestIter <- merge(performance, best_perf)$order
      performance$order <- NULL
    }
    
    
    ## Based on the optimality criterion, select the tuning parameter(s)
    bestTune <- performance[bestIter, paramNames, drop = FALSE]
  } else {
    bestTune <- tuneGrid
    performance <- evalSummaryFunction(y_dat, wts = weights,
                                       ctrl = trControl,
                                       lev = classLevels,
                                       metric = metric,
                                       method = method)
    perfNames <- names(performance)
    performance <- as.data.frame(t(performance), stringsAsFactors = TRUE)
    performance <- cbind(performance, tuneGrid)
    performance <- performance[-1,,drop = FALSE]
    tmp <- resampledCM <- NULL
    
  } # end(trControl$method != "none")
  
  ## Save some or all of the resampling summary metrics
  if(!(trControl$method %in% c("LOOCV", "oob", "none"))) {
    byResample <- switch(trControl$returnResamp,
                         none = NULL,
                         all = {
                           out <- resampleResults
                           colnames(out) <- gsub("^\\.", "", colnames(out))
                           out
                         },
                         final = {
                           out <- merge(bestTune, resampleResults)
                           out <- out[,!(names(out) %in% names(tuneGrid)), drop = FALSE]
                           out
                         })
  } else {
    byResample <- NULL
  }
  
  ## Reorder rows of performance
  orderList <- list()
  for(i in seq(along = paramNames)) orderList[[i]] <- performance[,paramNames[i]]
  
  performance <- performance[do.call("order", orderList),]
  
  if(trControl$verboseIter) {
    bestText <- paste(paste(names(bestTune), "=",
                            format(bestTune, digits = 3)),
                      collapse = ", ")
    if(nrow(performance) == 1) bestText <- "final model"
    cat("Fitting", bestText, "on full training set\n")
    flush.console()
  }
  
  ## Make the final model based on the tuning results
  indexFinal <- if(is.null(trControl$indexFinal))
    seq(along = data[[y_orig_val]]) else trControl$indexFinal
  
  if(!(length(trControl$seeds) == 1 && is.na(trControl$seeds)))
    set.seed(trControl$seeds[[length(trControl$seeds)]][1])
  finalTime <- system.time(
    finalModel <- rec_model(x,
                            subset_x(data, indexFinal),
                            method = models,
                            tuneValue = bestTune,
                            obsLevels = classLevels,
                            last = TRUE,
                            classProbs = trControl$classProbs,
                            sampling = trControl$sampling,
                            ...)
  )
  
  if(trControl$trim && !is.null(models$trim)) {
    if(trControl$verboseIter) old_size <- object.size(finalModel$fit)
    finalModel$fit <- models$trim(finalModel$fit)
    if(trControl$verboseIter) {
      new_size <- object.size(finalModel$fit)
      reduction <- format(old_size - new_size, units = "Mb")
      if(reduction == "0 Mb") reduction <- "< 0 Mb"
      p_reduction <- (unclass(old_size) - unclass(new_size))/unclass(old_size)*100
      p_reduction <- if(p_reduction < 1) "< 1%" else paste0(round(p_reduction, 0), "%")
      cat("Final model footprint reduced by", reduction, "or", p_reduction, "\n")
    }
  }
  
  trained_rec <- finalModel$recipe
  finalModel  <- finalModel$fit
  
  ## Remove this and check for other places it is reference
  ## replaced by tuneValue
  if(method == "pls") finalModel$bestIter <- bestTune
  
  ## To use predict.train and automatically use the optimal lambda,
  ## we need to save it
  if(method == "glmnet") finalModel$lambdaOpt <- bestTune$lambda
  
  if(trControl$returnData) {
    outData <- data
  } else outData <- NULL
  
  if(trControl$savePredictions == "final")
    tmp$predictions <- merge(bestTune, tmp$predictions)
  
  endTime <- proc.time()
  times <- list(everything = endTime - startTime,
                final = finalTime)
  
  out <- structure(list(method = method,
                        modelInfo = models,
                        modelType = modelType,
                        recipe = trained_rec,
                        results = performance,
                        pred = tmp$predictions,
                        bestTune = bestTune,
                        call = funcCall,
                        dots = list(...),
                        metric = metric,
                        control = trControl,
                        finalModel = finalModel,
                        trainingData = outData,
                        resample = byResample,
                        resampledCM = resampledCM,
                        perfNames = perfNames,
                        maximize = maximize,
                        yLimits = trControl$yLimits,
                        times = times,
                        levels = classLevels,
                        rs_seed = rs_seed),
                   class = c("train.recipe", "train"))
  trControl$yLimits <- NULL
  
  if(trControl$timingSamps > 0) {
    pData <- x_dat[sample(1:nrow(x_dat), trControl$timingSamps, replace = TRUE),,drop = FALSE]
    out$times$prediction <- system.time(predict(out, pData))
  } else  out$times$prediction <- rep(NA, 3)
  out
}


#' @method summary train
#' @export
summary.train <- function(object, ...) summary(object$finalModel, ...)

#' @importFrom stats predict residuals
#' @export
residuals.train <- function(object, ...) {
  if(object$modelType != "Regression") stop("train() only produces residuals on numeric outcomes", call. = FALSE)
  resid <- residuals(object$finalModel, ...)
  if(is.null(resid)) {
    if(!is.null(object$trainingData))  {
      resid <- object$trainingData$.outcome - predict(object, object$trainingData[, names(object$trainingData) != ".outcome",drop = FALSE])
    } else stop("The training data must be saved to produce residuals", call. = FALSE)
  }
  resid
}

#' @importFrom stats predict fitted
#' @export
fitted.train <- function(object, ...) {
  prd <- fitted(object$finalModel)
  if(is.null(prd)) {
    if(!is.null(object$trainingData)) {
      prd <- predict(object, object$trainingData[, names(object$trainingData) != ".outcome",drop = FALSE])
    } else stop("The training data must be saved to produce fitted values", call. = FALSE)
  }
  prd
  
}