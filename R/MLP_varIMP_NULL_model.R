MLP_varIMP_NULL_model=function(
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
  # Compute baseline metric for comparison
 # x_MLP=kerasR::expand_dims(train_x,axis = 2) 
  baseline <- as.data.frame(t(caret::postResample(
    pred = keras::predict_on_batch(optmodel, train_x),
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
  
  MLP_varIMP <- replicate(nsim,(plyr::llply(feature_names, .progress = "none",
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
                                              #x_tensor=kerasR::expand_dims(train_x_NULL,axis = 2) 
                                              permuted <- as.data.frame(t(caret::postResample(
                                                pred = keras::predict_on_batch(optmodel, train_x_NULL),
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
  #   "Importance" = apply(MLP_varIMP, MARGIN = 1, FUN = mean)
  #  )
  # if (nsim > 1) {
  #   tib$StDev <- apply(MLP_varIMP, MARGIN = 1, FUN = stats::sd)
  # }
  
  # Add all nsim scores as an attribute
  #  if (nsim > 1 && keep) {
  #   rownames(MLP_varIMP) <- feature_names
  #   colnames(MLP_varIMP) <- paste0("permutation_", seq_len(ncol(MLP_varIMP)))
  #   attr(tib, which = "raw_scores") <- MLP_varIMP
  #}
  
  # Add variable importance type attribute
  # attr(tib, which = "type") <- "permutation"
  
  # Add "vip" class
  # class(tib) <- c("MLP_IMP", class(tib))
  
  # Return results
  # tib
  varimp1= apply(MLP_varIMP, MARGIN = 1, function(x) do.call(cbind,x))
  varimp2=do.call(rbind,varimp1)
  ## now this is a array 1000, 30 dataframe
  nm <- colnames(varimp2)
  RMSE=varimp2[,grepl("^RMSE", nm)]
  Rsquared=varimp2[,grepl("^Rsquared", nm)]
  MAE=varimp2[,grepl("^MAE", nm)]
  
  RMSE=as.data.frame(rowMeans(RMSE,na.rm = TRUE))
  Rsquared=as.data.frame(rowMeans(Rsquared,na.rm = TRUE))
  MAE=as.data.frame(rowMeans(MAE,na.rm = TRUE))
  rownames(RMSE)=feature_names
  rownames(Rsquared)=feature_names
  rownames(MAE)=feature_names
  MLP_SNPsIMP=cbind(RMSE,Rsquared,MAE)
  rownames(MLP_SNPsIMP)=feature_names
  decrease_acc=lapply(1:3, function(i) baseline[,i]-MLP_SNPsIMP[,i])
  MLP_Decrease_acc=do.call(cbind,decrease_acc)
  return(list(MLP_Decrease_acc=MLP_Decrease_acc,MLP_SNPsIMP=MLP_SNPsIMP,baseline=baseline))
  
}

