### used in CNN and MLP with keras and Tensorflow package

"varImp.train" <- function(object, useModel = TRUE, nonpara = TRUE, scale = TRUE, ...) {
  code <- object$modelInfo
  if(is.null(code$varImp)) useModel <- FALSE
  if(useModel) {
    checkInstall(code$library)
    for(i in seq(along = code$library))
      do.call("requireNamespaceQuietStop", list(package = code$library[i]))
    imp <- code$varImp(object$finalModel, ...)
    modelName <- object$method
  } else {
    if(inherits(object, "train.recipe")) {
      x_dat <- recipes::juice(object$recipe, all_predictors())
      x_dat <- as.data.frame(x_dat, stringsAsFactors = FALSE)
      y_dat <- recipes::juice(object$recipe, all_outcomes())
      y_dat <- getElement(y_dat, names(y_dat))
    } else {
      isX <- which(!(colnames(object$trainingData) %in% ".outcome"))
      x_dat <- object$trainingData[, isX,drop = FALSE]
      y_dat <- object$trainingData[, -isX]
    }
    imp <- filterVarImp(x_dat, y_dat,
                        nonpara = nonpara,
                        ...)
    modelName <- ifelse(is.factor(y_dat),
                        "ROC curve",
                        ifelse(nonpara, "loess r-squared", "Linear model"))
  }
  if(scale) {
    if(class(object$finalModel)[1] == "pamrtrained") imp <- abs(imp)
    imp <- imp - min(imp, na.rm = TRUE)
    imp <- imp/max(imp, na.rm = TRUE)*100
  }
  out <- list(importance = imp,
              model = modelName,
              calledFrom = "varImp")
  
  structure(out, class = "varImp.train")
}



"varImp" <-
  function(object, ...){
    UseMethod("varImp")
  }
varImp.keras = function(object, ...) {
  object=keras::unserialize_model(object$finalModel$object)
  CNN_null_imp <- as.data.frame(CNN_varIMP_NULL_model(object),
                                             feature_names,
                                             train_y=env,
                                             train_x=as.matrix(genotype),
                                             #smaller_is_better = FALSE,
                                             type = c("difference", "ratio"),
                                             nsim = 100,## MCMC permutation large numbers need much more time
                                             sample_size = NULL,
                                             sample_frac = NULL,
                                             verbose = FALSE,
                                             progress = "none",
                                             parallel = TRUE,
                                             paropts = NULL)
  out=CNN_null_imp$CNN_Decrease_acc
  colnames(out)[colnames(out) == "CNN_Decrease_acc"] <- "Overall"
  rownames(out) <- out$variable
  out[, c("Overall"), drop = FALSE]
  
}




#' @importFrom stats4 coef
GarsonWeights <- function(object)
{
  beta <- coef(object)
  abeta <- abs(beta)
  nms <- names(beta)
  i2h <- array(NA, dim = object$n[2:1])
  h2o <- array(NA, dim = object$n[2:3])
  
  for (hidden in 1:object$n[2]) {
    for (input in 1:object$n[1]) {
      label <- paste("i", input, "->h", hidden,"$", sep = "")
      i2h[hidden, input] <- abeta[grep(label, nms, fixed = FALSE)]
    }
  }
  for(hidden in 1:object$n[2]){
    for(output in 1:object$n[3]){
      label <- paste("h", hidden, "->o",
                     ifelse(object$n[3] == 1, "", output),
                     sep = "")
      h2o[hidden,output] <- abeta[grep(label, nms, fixed = TRUE)]
    }
  }
  
  if(FALSE)
  {
    ## Test case from Gevrey, M., Dimopoulos, I., & Lek,
    ## S. (2003). Review and comparison of methods to study the
    ## contribution of variables in artificial neural network
    ## models. ecological modelling, 160(3), 249-264.
    i2h <- matrix(c(-1.67624,  3.29022,  1.32466,
                    -0.51874, -0.22921, -0.25526,
                    -4.01764,  2.12486, -0.08168,
                    -1.75691, -1.44702,  0.58286),
                  ncol = 3, byrow = TRUE)
    h2o <- matrix(c(4.57857, -0.48815, -5.73901, -2.65221),
                  ncol = 1)
  }
  
  ##  From Gevrey et al. (2003): "For each hidden neuron i, multiply
  ##  the absolute value of the hidden-output layer connection
  ##  weight by the absolute value of the hidden-input layer
  ##  connection weight. Do this for each input variable j. The
  ##  following products Pij are obtained"
  
  
  ## We'll do this one response at a time. Gevrey et al. (2003) do
  ## not discuss multiple outputs, but the results are the same (at
  ## least in the case of classification).
  
  imp <- matrix(NA, nrow = object$n[1], ncol = object$n[3])
  
  
  for(output in 1:object$n[3])
  {
    Pij <- i2h * NA
    for(hidden in 1:object$n[2]) Pij[hidden,] <- i2h[hidden,] * h2o[hidden,output]
    
    ## "For each hidden neuron, divide Pij by the sum for all the
    ## input variables to obtain Qij. For example for Hidden 1, Q11 =
    ## P11/(P11+P12+P13).
    
    Qij <- Pij * NA
    for(hidden in 1:object$n[2]) Qij[hidden,] <- Pij[hidden,] / sum(Pij[hidden,])
    
    
    ## "For each input neuron, sum the product Sj formed from the
    ## previous computations of Qij. For example, S1 =
    ## Q11+Q21+Q31+Q41."
    
    Sj <- apply(Qij, 2, sum)
    
    ## "Divide Sj by the sum for all the input variables. Expressed as
    ## a percentage, this gives the relative importance or
    ## distribution of all output weights attributable to the given
    ## input variable. For example, for the input neuron 1, the
    ## relative importance is equal to (S1/100)/(S1+S2+S3)"
    
    imp[,output] <- Sj/sum(Sj)*100
    rm(Pij, Qij, Sj)
  }
  
  colnames(imp) <- if(!is.null(colnames(object$residuals))) colnames(object$residuals) else paste("Y", 1:object$n[3], sep = "")
  rownames(imp) <- if(!is.null(object$coefnames)) object$coefnames else  paste("X", 1:object$n[1], sep = "")
  imp
}


GarsonWeights_FCNN4R <- function (object, xnames = NULL, ynames = NULL) {
  beta <- abs(object$net@m_w_values[which(object$net@m_w_flags != 0L)])
  dims <- object$net@m_layers
  
  index <- (dims[1]+1)*dims[2]
  i2h <- t(matrix(beta[1:index], ncol = dims[2]))
  i2h <- i2h[, -1,drop = FALSE]
  
  h2o <- matrix(beta[(index+1):length(beta)], ncol = dims[3])
  h2o <- h2o[-1,,drop = FALSE]
  
  imp <- matrix(NA, nrow = dims[1], ncol = dims[3])
  for (output in 1:dims[3]) {
    Pij <- i2h * NA
    for (hidden in 1:dims[2]) Pij[hidden, ] <- i2h[hidden,] * h2o[hidden, output]
    Qij <- Pij * NA
    for (hidden in 1:dims[2]) Qij[hidden, ] <- Pij[hidden,]/sum(Pij[hidden, ])
    Sj <- apply(Qij, 2, sum)
    imp[, output] <- Sj/sum(Sj) * 100
    rm(Pij, Qij, Sj)
  }
  rownames(imp) <- if(is.null(xnames))
    paste("X", 1:dims[1], sep = "") else
      xnames
  colnames(imp) <- if(is.null(ynames))
    paste("Y", 1:dims[3], sep = "") else
      ynames
  imp
}

varImpDependencies <- function(libName){
  code <- getModelInfo(libName, regex = FALSE)[[1]]
  checkInstall(code$library)
  for(i in seq(along = code$library))
    do.call("requireNamespaceQuietStop", list(package = code$library[i]))
  return(code)
}





varImp.neuralnet <- function(object, ...){
  code <- varImpDependencies("neuralnet")
  code$varImp(object, ...)
}
varImp.RSNNS <- function(object, ...){
  code <- varImpDependencies("RSNNS")
  code$varImp(object, ...)
}
varImp.h2o <- function(object, ...){
  code <- varImpDependencies("h2o")
  code$varImp(object, ...)
}
varImp.keras <- function(object, ...){
  code <- varImpDependencies("kears")
  code$varImp(object, ...)
}
#' @rdname varImp
#' @export
varImp.bagEarth <- function(object, ...){
  code <- varImpDependencies("bagEarth")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.bagFDA <- function(object, ...){
  code <- varImpDependencies("bagFDA")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.C5.0 <- function(object, ...){
  code <- varImpDependencies("C5.0")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.cubist <- function(object, weights = c(0.5, 0.5), ...){
  code <- varImpDependencies("cubist")
  code$varImp(object, weights = weights, ...)
}

#' @rdname varImp
#' @export
varImp.dsa <- function(object, cuts = NULL, ...){
  code <- varImpDependencies("partDSA")
  code$varImp(object, cuts = cuts, ...)
}

#' @rdname varImp
#' @export
varImp.glm <- function(object, ...){
  code <- varImpDependencies("glm")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.glmnet <- function(object, lambda = NULL, ...){
  code <- varImpDependencies("glmnet")
  code$varImp(object, lambda = lambda, ...)
}

#' @rdname varImp
#' @export
varImp.JRip <- function(object, ...){
  code <- varImpDependencies("JRip")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.multinom <- function(object, ...){
  code <- varImpDependencies("multinom")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.nnet <- function(object, ...){
  code <- varImpDependencies("nnet")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.avNNet <- function(object, ...){
  code <- varImpDependencies("nnet")
  imps <- lapply(object$model, code$varImp)
  imps <- do.call("rbind", imps)
  imps <- aggregate(imps, by = list(vars = rownames(imps)), mean)
  rownames(imps) <- as.character(imps$vars)
  imps$vars <- NULL
  imps
}

foo <- function(x) {
  browser()
  colMeans(x)
}


#' @rdname varImp
#' @export
varImp.PART <- function(object, ...){
  code <- varImpDependencies("PART")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.RRF <- function(object, ...){
  code <- varImpDependencies("RRF")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.rpart <- function(object, surrogates = FALSE, competes = TRUE, ...){
  code <- varImpDependencies("rpart")
  code$varImp(object, surrogates = surrogates, competes = competes, ...)
}

#' @rdname varImp
#' @export
varImp.randomForest <- function(object, ...){
  code <- varImpDependencies("rf")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.gbm <- function(object, numTrees = NULL, ...){
  code <- varImpDependencies("gbm")
  code$varImp(object, numTrees = numTrees, ...)
}

#' @rdname varImp
#' @export
varImp.classbagg <- function(object, ...){
  code <- varImpDependencies("treebag")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.regbagg <- function(object, ...){
  code <- varImpDependencies("treebag")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.pamrtrained <- function(object, threshold, data, ...){
  code <- varImpDependencies("pam")
  code$varImp(object,
              threshold = object$bestTune$threshold,
              data = object$finalModel$xData,
              ...)
}

#' @rdname varImp
#' @export
varImp.lm <- function(object, ...){
  code <- varImpDependencies("lm")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.mvr <- function(object, estimate = NULL, ...){
  code <- varImpDependencies("pls")
  code$varImp(object, estimate = estimate, ...)
}

#' @rdname varImp
#' @export
varImp.earth <- function(object, value = "gcv", ...){
  code <- varImpDependencies("earth")
  code$varImp(object, value = value, ...)
}

#' @rdname varImp
#' @export
varImp.RandomForest <- function(object, ...){
  code <- varImpDependencies("cforest")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.plsda <- function(object, ...){
  code <- varImpDependencies("pls")
  code$varImp(object, ...)
}

#' @rdname varImp
#' @export
varImp.fda <- function(object, value = "gcv", ...){
  code <- varImpDependencies("fda")
  code$varImp(object, value = value, ...)
}

#' @rdname varImp
#' @export
varImp.gam <- function(object, ...){
  code <- varImpDependencies("gam")
  code$varImp(object, ...)
}


#' @rdname varImp
#' @export
varImp.Gam <- function(object, ...){
  code <- varImpDependencies("gamSpline")
  code$varImp(object, ...)
}

