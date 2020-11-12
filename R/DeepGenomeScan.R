#####DeepGenomeScan 

"DeepGenomeScan" <-
    function(genotype, ...){
      UseMethod("DeepGenomeScan")
    }


DeepGenomeScan.default=function(genotype, env,method = "mlp",
                                metric="RMSE",
                                preProcess = NULL,
                                trControl =trainControl(## 5-fold CV, repeat 5 times
                                  method = "repeatedcv",
                                  number = 5,
                                  ## repeated 5 times
                                  repeats = 5,search = "random"),
                                 tuneLength = 10,...){ 
  model_train=caret::train(x=genotype,
                           y=env,
                           method = method,
                           preProcess = preProcess,
                           metric =metric ,
                           trControl = trControl,
                           tuneLength = tuneLength,...)
  class(model_train)="DeepGenomeScan"
  return(model_train)  
  
  }

DeepGenomeScan.CNN=function(genotype, env,method = "mlp",
                        metric = "RMSE",
                         preProcess = NULL,
                        trControl = trainControl(## 5-fold CV, repeat 5 times
                          method = "repeatedcv",
                          number = 5,
                          ## repeated 5 times
                          repeats = 5,search = "random"),
                        tuneLength = 10,
                        importance = c("permutation","NULL_model","Olden","Garson"), nsim= if (importance=="permutation") nsim=10,...){ 
         
    model_train=caret::train(x=genotype,
                                y=env,
                             method = method,
                             preProcess = preProcess,
                             metric =metric ,
                             trControl = trControl,
                             tuneLength = tuneLength, ... )
  
  
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
  class(model_train)="DeepGenomeScan"
  return(model_train)
}

DeepGenomeScan.formula=function(form, data,...){

  model_train=caret::train(form=form, data=data, ...)
  class(model_train)="DeepGenomeScan"
  return(model_train)
}

DeepGenomeScan.recipe=function(genotype, data, method = "mlp", metric = "RMSE",
                                 
                                 trControl = trainControl(),
                                 
                                 tuneLength = ifelse(trControl$method == "none", 1, 3),...){
  model_train=caret::train(x=genotype,data=data, method = method,
                            metric =metric ,
                           trControl = trControl,
                           tuneLength = tuneLength, ...)
  class(model_train)="DeepGenomeScan"
  return(model_train)
}

#### Improtance function





DeepGenomeScan.keras=function(genotype, env,method=method,
                        metric = "MAE",## "Accuracy", "RMSE","Rsquared"
                        preProcess = NULL,
                          tuneLength = 11, ### number of parameter combinations
                        # tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                        trControl = trainControl(## 5-fold CV, repeat 5 times
                          method = "repeatedcv",
                          number = 5,
                          ## repeated 5 times
                          repeats = 5,search = "random"),
                        importance = c("permutation","NULL_model","Olden","Garson"), nsim= if (importance=="permutation") nsim=10,...){
  

 DLmodel<- caret::train(x=genotype,y=env,
                                 method=method,
                                 metric = metric,## "Accuracy", "RMSE","Rsquared"
                        preProcess = preProcess,
                                 tuneLength = tuneLength, ### 11 tunable parameters 11^2
                                 # tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                                 ...,
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
#load("sim_example.RData")
#genotype=sim_example[,-c(1:14)]
#env=sim_example[,2:11]
