#library(vita)

### make sure the plyr
#nsim=10 ####the number of Monte Carlo replications to perform. Default is 1. If nsim > 1, the results from each replication are simply averaged together (the standard deviation will also be returned).
#feature_names=colnames(data1[,-1])
#parallel = FALSE
#paropts = NULL
#sample_size=NULL
#train_y=as.matrix(data1[,1])
#train_x= data1[,-1]
#metric = NULL
#optmodel=keras::unserialize_model(model_Keras_CNN$finalModel$object)

#set.seed(123)
#CNN_sim_imp=CNN_varIMP_permute(optmodel,
#                               feature_names,
 #                              train_y,
#                               train_x,
#                               #smaller_is_better = FALSE,
#                               type = c("difference", "ratio"),
#                               nsim = 1,## large numbers need much more time
#                               sample_size = NULL,
#                               sample_frac = NULL,
#                               verbose = FALSE,
 #                              progress = "none",
#                               parallel = TRUE,
#                               paropts = NULL)

#CNN_sim_impNULL=CNN_varIMP_NULL_model(optmodel,
#                                      feature_names,
#                                      train_y,
 #                                     train_x,
                                      #smaller_is_better = FALSE,
#                                      type = c("difference", "ratio"),
#                                      nsim = 1,## large numbers need much more time
 #                                     sample_size = NULL,
#                                      sample_frac = NULL,
##                                      verbose = FALSE,
 #                                     progress = "none",
#                                      parallel = TRUE,
#                                      paropts = NULL)
#Sys.time()
#cor.test(CNN_sim_impNULL$CNN_Decrease_acc[,1],CNN_sim_imp$CNN_Decrease_acc[,1])

CNN_varIMP_permute <- function(
  optmodel,
  feature_names = NULL,
  train_y = NULL,
  train_x = NULL,
  #smaller_is_better = NULL,
  type = c("difference", "ratio"),
  nsim = 1,
  sample_size = NULL,
  sample_frac = NULL,
  verbose = FALSE,
  progress = "none",
  parallel = FALSE,
  paropts = NULL,
  ...
) {
  require(plyr)
  require(kerasR)
  require(keras)
  # Compute baseline metric for comparison
  x_cnn=kerasR::expand_dims(train_x,axis = 2) 
  baseline <- as.data.frame(t(caret::postResample(
    pred = keras::predict_on_batch(optmodel, x_cnn),
    obs = train_y)))
  (rm(x_cnn))
  # Type of comparison
  type <- match.arg(type)
  `%compare%` <- if (type == "difference") {
    `-`
  } else {
    `/`
  }
  
  # Construct VI scores
  #
  # Loop through each feature and do the following:
  #
  #   1. make a copy of the training data;
  #   2. permute the values of the original feature;
  #   3. get new predictions based on permuted data set;
  #   4. record difference in accuracy.
  
  
  permute_columns <- function(x, columns = NULL) {
    if (is.null(columns)) {
      stop("No columns specified for permutation.")
    }
    x[, columns] <- x[sample(nrow(x)), columns]
    x
  }
  
  
  #
  sort_importance_scores <- function(x, decreasing) {
    x[order(x$Importance, decreasing = decreasing), ]
  }
  #unlist
  CNN_varIMP <- replicate(nsim,(plyr::llply(feature_names, .progress = "none",
                                            .parallel = parallel, .paropts = paropts,
                                            .fun = function(x) {
                                              if (verbose && !parallel) {
                                                message("Computing variable importance for ", x, "...")
                                              }
                                              if (!is.null(sample_size)) {
                                                ids <- sample(length(train_y), size = sample_size, replace = FALSE)
                                                train_x <- train_x[ids, ]
                                                train_y <- train_y[ids]
                                              }
                                              # train_x_permuted <- train_x  # make copy
                                              # train_x_permuted[[x]] <- sample(train_x_permuted[[x]])  # permute values
                                              train_x_permuted <- permute_columns(train_x, columns = x)
                                              ###caret::postResample(pred, data1[,1]) give three values: RMSE Rsquared, MAE
                                              x_tensor=kerasR::expand_dims(train_x_permuted,axis = 2) 
                                              permuted <- as.data.frame(t(caret::postResample(
                                                pred = keras::predict_on_batch(optmodel, x_tensor),
                                                obs = train_y)))
                                              # if (smaller_is_better) {
                                              #   permuted %compare% baseline  # e.g., RMSE
                                              # } else {
                                              #    baseline %compare% permuted  # e.g., R-squared
                                              # }
                                            })
  ))
  
  # Construct tibble of variable importance scores
  #tib <- tibble::tibble(
  #  "Variable" = feature_names,
  #  "Importance" = apply(CNN_varIMP, MARGIN = 1, FUN = mean)
  #)
  #if (nsim > 1) {
  #  tib$StDev <- apply(CNN_varIMP, MARGIN = 1, FUN = stats::sd)
  #}
  
  # Add all nsim scores as an attribute
  #if (nsim > 1 && keep) {
  # rownames(CNN_varIMP) <- feature_names
  #  colnames(CNN_varIMP) <- paste0("permutation_", seq_len(ncol(CNN_varIMP)))
  #  attr(tib, which = "raw_scores") <- CNN_varIMP
  #}
  
  # Add variable importance type attribute
  #attr(tib, which = "type") <- "permutation"
  
  # Add "vip" class
  #class(tib) <- c("CNN_IMP", class(tib))
  
  # Return results
  #tib#
  CNN_SNPsIMP=do.call(rbind,CNN_varIMP)
  rownames(CNN_SNPsIMP)=feature_names
  decrease_acc=lapply(1:3, function(i) baseline[,i]-CNN_SNPsIMP[,i])
  CNN_Decrease_acc=do.call(cbind,decrease_acc)
  return(list(CNN_Decrease_acc=CNN_Decrease_acc,CNN_SNPsIMP=CNN_SNPsIMP,baseline=baseline))
  
}


###########################one risk is the purmution method is not stable and the CNN pooled layer may decode the complex relations. Therefore, the 
CNN_varIMP_NULL_model <- function(
  optmodel,
  feature_names = NULL,
  train_y = NULL,
  train_x = NULL,
  smaller_is_better = NULL,
  type = c("difference", "ratio"),
  nsim = 1,
  sample_size = NULL,
  sample_frac = NULL,
  verbose = FALSE,
  progress = "none",
  parallel = FALSE,
  paropts = NULL,
  ...
) {
  require(caret)
  require(plyr)
  require(kerasR)
  require(keras)
  # Compute baseline metric for comparison
  x_cnn=kerasR::expand_dims(train_x,axis = 2) 
  baseline <- as.data.frame(t(caret::postResample(
    pred = keras::predict_on_batch(optmodel, x_cnn),
    obs = train_y)))
  
  # Type of comparison
  type <- match.arg(type)
  `%compare%` <- if (type == "difference") {
    `-`
  } else {
    `/`
  }
  
  # Construct VIP scores
  #
  # Loop through each feature and do the following:
  #
  #   1. make a copy of the training data;
  #   2. permute the values of the original feature;
  #   3. get new predictions based on permuted data set;
  #   4. record difference in accuracy.
  
  
  NULL_columns <- function(x, columns = NULL) {
    if (is.null(columns)) {
      stop("No columns specified for permutation.")
    }
    x[, columns] <- rep(1,nrow(x))
    x
  }
  
  
  #
  sort_importance_scores <- function(x, decreasing) {
    x[order(x$Importance, decreasing = decreasing), ]
  }
  
  CNN_varIMP <- replicate(nsim,(plyr::llply(feature_names, .progress = "none",
                                            .parallel = parallel, .paropts = paropts,
                                            .fun = function(x) {
                                              if (verbose && !parallel) {
                                                message("Computing variable importance for ", x, "...")
                                              }
                                              if (!is.null(sample_size)) {
                                                ids <- sample(length(train_y), size = sample_size, replace = FALSE)
                                                train_x <- train_x[ids, ]
                                                train_y <- train_y[ids]
                                              }
                                              # train_x_permuted <- train_x  # make copy
                                              # train_x_permuted[[x]] <- sample(train_x_permuted[[x]])  # permute values
                                              train_x_NULL <- NULL_columns(train_x, columns = x)
                                              ###caret::postResample(pred, data1[,1]) give three values: RMSE Rsquared, MAE
                                              x_tensor=kerasR::expand_dims(train_x_NULL,axis = 2) 
                                              permuted <- as.data.frame(t(caret::postResample(
                                                pred = keras::predict_on_batch(optmodel, x_tensor),
                                                obs = train_y)))
                                              #  if (smaller_is_better) {
                                              #     permuted %compare% baseline  # e.g., RMSE
                                              #  } else {
                                              #    baseline %compare% permuted  # e.g., R-squared
                                              #   }
                                            })
  ))
  
  # Construct tibble of variable importance scores
  #  tib <- tibble::tibble(
  #    "Variable" = feature_names,
  #   "Importance" = apply(CNN_varIMP, MARGIN = 1, FUN = mean)
  #  )
  # if (nsim > 1) {
  #   tib$StDev <- apply(CNN_varIMP, MARGIN = 1, FUN = stats::sd)
  # }
  
  # Add all nsim scores as an attribute
  #  if (nsim > 1 && keep) {
  #   rownames(CNN_varIMP) <- feature_names
  #   colnames(CNN_varIMP) <- paste0("permutation_", seq_len(ncol(CNN_varIMP)))
  #   attr(tib, which = "raw_scores") <- CNN_varIMP
  #}
  
  # Add variable importance type attribute
  # attr(tib, which = "type") <- "permutation"
  
  # Add "vip" class
  # class(tib) <- c("CNN_IMP", class(tib))
  
  # Return results
  # tib
  CNN_SNPsIMP=do.call(rbind,CNN_varIMP)
  rownames(CNN_SNPsIMP)=feature_names
  decrease_acc=lapply(1:3, function(i) baseline[,i]-CNN_SNPsIMP[,i])
  CNN_Decrease_acc=do.call(cbind,decrease_acc)
  return(list(CNN_Decrease_acc=CNN_Decrease_acc,CNN_SNPsIMP=CNN_SNPsIMP,baseline=baseline))
  
}




