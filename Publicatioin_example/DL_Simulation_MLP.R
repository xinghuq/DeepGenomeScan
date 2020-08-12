library(caret)### for ML calling functions and performance estimation
library(keras) ### for DL
library("tensorflow")
library("caretEnsemble")
library(rBayesianOptimization) ### for parameter tuning and optimazation
#library("RSNNS")
set.seed(999)
options(warn=-1)

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

dirs=list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
fls <- list.files(dirs, pattern="*.csv", full.names = TRUE)
simdata=lapply(fls, read.csv)
#save(simdata,file = "Dimulated_Data.RData")
#### separate data, while this is not necrssary as long as you know and subset the data
env=normalize(simdata[[1]][,2:11])
genotype=simdata[[1]][,-(1:14)]
qtrait=simdata[[1]][,12:14]
para=colnames(env)

#### before run  the model, set parallel (only in the cluster, for my HP desktop, only 8 GB is not necessary to do this)
library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
### some cluster may not able to use the parallel when calling keras Tensorflow
### 
Sys.time()





########## MLP random search############################################
set.seed(123)
econtrol <- trainControl(## 5-fold CV, repeat 5 times
  method = "repeatedcv",
  number = 5,
  ## repeated ten times
  repeats = 5,search = "random")

print("Start traing DL mlp")
Sys.time()
for(j in 1:length(simdata)) {
  for (i in 1:length(para)){
    
    data0=(cbind(simdata[[j]][,para[i]],simdata[[j]][,-(1:14)]))
    # colnames(data0)[1]=para[[i]]
    colnames(data0)[1]="enviri" ### This is more efficient than calling the name of the env factor when traing the model
    data1=as.data.frame(apply(data0,2,normalize))
    model_mlpKerasDropout_envi<- train(enviri~., data=data1,
                                       method=c("mlpKerasDropout"),
                                       metric = "MAE",## "Accuracy", "RMSE"
                                       tuneLength = 100, ### search 100 combinations of parameters
                                       verbose=0,# verbose=1 is reporting the progress,o is sclience
                                       trControl = econtrol,importance = TRUE)
    
    #  keras::unserialize_model(gbmGrid_model_mlpKerasDropout_env$finalModel$object)
    MLPImp <- varImp(model_mlpKerasDropout_envi, scale = FALSE)
    write.csv(MLPImp$importance,file = paste0("sim_",j,"_",para[i],"_mlpKerasDropout_env.csv"))
    save(model_mlpKerasDropout_envi,file=paste0("sim_",j,"_",para[i],"_mlpKerasDropout_env.RData"))
    write.csv(model_mlpKerasDropout_envi$results,file = paste0("sim_",j,"_",para[i],"_model_mlpKerasDropout_envi_tuning.csv") )
    #registerDoSEQ()
    #stopCluster(cl)
  }
}



Sys.time()



save.image(file = "MLP_final_after_trained_data.RData")
