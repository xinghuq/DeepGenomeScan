#####DeepGenomeScan 

"DeepGenomeScan" <-
    function(x, ...){
      UseMethod("DeepGenomeScan")
    }

DeepGenomeScan.default=function(genotype, env,method = "mlpWeightDecayML",preProcess = NULL,weights = NULL,
                        metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                        maximize = ifelse(metric %in% c("RMSE", "logLoss", "MAE"), FALSE, TRUE),
                        trControl = trainControl(),
                        tuneGrid = NULL, tuneLength = ifelse(trControl$method == "none", 1, 3),
                        importance = c("permutation","NULL_model","Olden","Garson"), nsim= if (importance=="permutation") nsim=10,...){ 
  model_train=caret::train(x=as.matrix(genotype),y=as.matrix(env),method = method,preProcess = preProcess,..., weights = weight,metric =metric ,maximize =maximize,
    trControl = trControl,
    tuneGrid = tuneGrid,
    tuneLength = tuneLength)
  
  
  if (method!= "CNN_sgd" | "*Keras" | importance =="permutation" |"NULL_model")
    #  keras::unserialize_model(gbmGrid_model_mlpKerasDropout_env$finalModel$object)
    optmodel=keras::unserialize_model(model_train$finalModel$object)
  set.seed(123)
  feature_names=colnames(genotype)
  Tensor_sim_imp=CNN_varIMP_permute(optmodel,
                                    feature_names,
                                    train_y=env,
                                    train_x=genotype,
                                    #smaller_is_better = FALSE,
                                    type = c("difference", "ratio"),
                                    nsim = 10,## large numbers need much more time
                                    sample_size = NULL,
                                    sample_frac = NULL,
                                    verbose = FALSE,
                                    progress = "none",
                                    parallel = TRUE,
                                    paropts = NULL)
  VarImps=Tensor_sim_imp$CNN_Decrease_acc
  # write.csv(Tensor_sim_imp$CNN_Decrease_acc,file = "Tensor_sim_imp.csv")
  if (method!= "CNN_sgd" | "*Keras" | importance== "NULL_model")
    Tensor_sim_impNULL=CNN_varIMP_NULL_model(optmodel,
                                             feature_names,
                                             train_y=env,
                                             train_x=genotype,
                                             #smaller_is_better = FALSE,
                                             type = c("difference", "ratio"),
                                             nsim = 10,## large numbers need much more time
                                             sample_size = NULL,
                                             sample_frac = NULL,
                                             verbose = FALSE,
                                             progress = "none",
                                             parallel = TRUE,
                                             paropts = NULL)
  VarImps=Tensor_sim_impNULL$CNN_Decrease_acc
  #write.csv(Tensor_sim_impNULL$CNN_Decrease_acc,file = "Tensor_sim_impNULL.csv")
  if (importance=="Olden") VarImps <- NeuralNetTools::olden(model_train,bar_plot = FALSE)
  if (importance=="Garson") VarImps <- NeuralNetTools::garson(model_train,bar_plot = FALSE)
  
  return(list(model_train,VarImps))
  
  return(model_train)
}

DeepGenomeScan.formula=function(form, data,  weights = NULL, subset = NULL, na.action = na.fail, contrasts = NULL,...){

  model_train=caret::train(form, data, ..., weights, subset, na.action = na.action, contrasts = contrasts)
return(model_train)
}

DeepGenomeScan.recipe=function(x, data, method = "mlpWeightDecayML", metric = ifelse(is.factor(y_dat), "Accuracy", "RMSE"),
                                 maximize = ifelse(metric %in% c("RMSE", "logLoss", "MAE"), FALSE, TRUE),
                                 trControl = trainControl(),
                                 tuneGrid = NULL,
                                 tuneLength = ifelse(trControl$method == "none", 1, 3),...){
  
  model_train=caret::train(x,data, method = method,..., metric = metric,
    maximize = maximize,
    trControl = trControl,
    tuneGrid = tuneGrid,
    tuneLength = tuneLength)
  return(model_train)
}

#### Improtance function

NeuralNetTools::olden




DeepGenomeScan=function(genotype, env,method=modelmlpkerasdropout,
                        metric = "MAE",## "Accuracy", "RMSE","Rsquared"
                        tuneLength = 11, ### number of parameter combinations
                        # tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                        verbose=0,# verbose=1 is reporting the progress,o is sclience
                        trControl = trainControl(## 5-fold CV, repeat 5 times
                          method = "repeatedcv",
                          number = 5,
                          ## repeated 5 times
                          repeats = 5,search = "random"),
                        importance = c("permutation","NULL_model","Olden","Garson"), nsim= if (importance=="permutation") nsim=10,...){
  

 DLmodel<- caret::train(x=as.matrix(genotype),y=as.matrix(env),
                                 method=method,
                                 metric = metric,## "Accuracy", "RMSE","Rsquared"
                                 tuneLength = tuneLength, ### 11 tunable parameters 11^2
                                 # tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                                 verbose=verbose,# verbose=1 is reporting the progress,o is sclience
                                 trControl = trControl,importance = importance)
  
  if (method!= "CNN_sgd" | "*Keras" | importance =="permutation" |"NULL_model")
 #  keras::unserialize_model(gbmGrid_model_mlpKerasDropout_env$finalModel$object)
  optmodel=keras::unserialize_model(DLmodel$finalModel$object)
  set.seed(123)
  feature_names=colnames(genotype)
  Tensor_sim_imp=CNN_varIMP_permute(optmodel,
                                 feature_names,
                                 train_y=env,
                                 train_x=genotype,
                                 #smaller_is_better = FALSE,
                                 type = c("difference", "ratio"),
                                 nsim = 10,## large numbers need much more time
                                 sample_size = NULL,
                                 sample_frac = NULL,
                                 verbose = FALSE,
                                 progress = "none",
                                 parallel = TRUE,
                                 paropts = NULL)
  VarImps=Tensor_sim_imp$CNN_Decrease_acc
 # write.csv(Tensor_sim_imp$CNN_Decrease_acc,file = "Tensor_sim_imp.csv")
  if (method!= "CNN_sgd" | "*Keras" | importance== "NULL_model")
   Tensor_sim_impNULL=CNN_varIMP_NULL_model(optmodel,
                                        feature_names,
                                        train_y=env,
                                        train_x=genotype,
                                        #smaller_is_better = FALSE,
                                        type = c("difference", "ratio"),
                                        nsim = 10,## large numbers need much more time
                                        sample_size = NULL,
                                        sample_frac = NULL,
                                        verbose = FALSE,
                                        progress = "none",
                                        parallel = TRUE,
                                        paropts = NULL)
  VarImps=Tensor_sim_impNULL$CNN_Decrease_acc
  #write.csv(Tensor_sim_impNULL$CNN_Decrease_acc,file = "Tensor_sim_impNULL.csv")
   if (importance=="Olden") VarImps <- NeuralNetTools::olden(DLmodel.bar_plot = FALSE)
   if (importance=="Garson") VarImps <- NeuralNetTools::garson(DLmodel.bar_plot = FALSE)
  
return(list(DLmodel,VarImps))
  
}


#### example
load("sim_example.RData")
genotype=sim_example[,-c(1:14)]
env=sim_example[,2:11]
