library(caret)### for ML calling functions and performance estimation
library(keras) ### for DL
library("tensorflow")
library("caretEnsemble")
library(kerasR)
#library("RSNNS")
set.seed(123)
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
#library(doParallel)
#cl <- makePSOCKcluster(5)
#registerDoParallel(cl)
### some cluster may not able to use the parallel when calling keras Tensorflow
### 
Sys.time()

########## CNN random search############################################
set.seed(123)
econtrol <- trainControl(## 5-fold CV, repeat 5 times
  method = "repeatedcv",
  number = 5,
  ## repeated ten times
  repeats = 5,search = "random")
### parameters for CNN importance
nsim=1####the number of Monte Carlo replications to perform. Default is 1. If nsim > 1, the results from each replication are simply averaged together (the standard deviation will also be returned).
feature_names=colnames(genotype)
print("Start traing DL CNN")
Sys.time()
for(j in 1:length(simdata)) {
  for (i in 1:length(para)){
    print(paste0("sim_",j,"_",para[i],"tuning CNN 1,000 parameters"))
    Sys.time()
    data0=(cbind(simdata[[j]][,para[i]],simdata[[j]][,-(1:14)]))
    # colnames(data0)[1]=para[[i]]
    colnames(data0)[1]="enviri" ### This is more efficient than calling the name of the env factor when traing the model
    #data1=as.data.frame(apply(data0,2,normalize))
    model_CNN_L_Dropout_envi<- train(enviri~., data=data0,
                                       method=CNN_sgd,
                                       metric = "MAE",## "Accuracy", "RMSE"
                                       preProcess=c("scale","center"),
                                       tuneLength = 1000, ### search 500 combinations of parameters
                                       verbose=0,# verbose=1 is reporting the progress,o is sclience
                                       trControl = econtrol,importance = TRUE)
    print(paste0("sim_",j,"_",para[i],"tuning CNN finished"))
    Sys.time()
    #  keras::unserialize_model(gbmGrid_model_CNN_L_Dropout_env$finalModel$object)
    ### loess regression importance
    CNNImploess <- varImp(model_CNN_L_Dropout_envi, scale = FALSE)
    write.csv(CNNImploess$importance,file = paste0("sim_",j,"_",para[i],"_CNN_L_Dropout_loess_importance_env.csv"))
    save(model_CNN_L_Dropout_envi,file=paste0("sim_",j,"_",para[i],"_CNN_L_Dropout_env.RData"))
    write.csv(model_CNN_L_Dropout_envi$results,file = paste0("sim_",j,"_",para[i],"_model_CNN_L_Dropout_envi_tuning.csv") )
    ###CNN importance 
    Sys.time()
    print(paste0("sim_",j,"_",para[i],"CNN optimal model permuting importance"))
    optmodel=keras::unserialize_model(model_CNN_L_Dropout_envi$finalModel$object)
    ### make sure the plyr
    set.seed(123)
    train_y=as.matrix(data0$enviri)
    train_x= data0[,-1]
    CNN_sim_imppermute=CNN_varIMP_permute(optmodel,
                                   feature_names,
                                   train_y,
                                   train_x,
                                   #smaller_is_better = FALSE,
                                   type = c("difference", "ratio"),
                                   nsim = 1,## large numbers need much more time
                                   sample_size = NULL,
                                   sample_frac = NULL,
                                   verbose = FALSE,
                                   progress = "none",
                                   parallel = TRUE,
                                   paropts = NULL)
    save(CNN_sim_imppermute,file = paste0("sim_",j,"_",para[i],"_model_CNN_L_Dropout_imppermute.RData") )
    write.csv(CNN_sim_imppermute$CNN_Decrease_acc,file =paste0("sim_",j,"_",para[i],"_model_CNN_L_Dropout_imppermute_Decrease_acc.csv") )
    Sys.time()
    print(paste0("sim_",j,"_",para[i],"CNN optimal model NULL model permuting importance"))
    CNN_sim_impNULL=CNN_varIMP_NULL_model(optmodel,
                                          feature_names,
                                          train_y,
                                          train_x,
                                          #smaller_is_better = FALSE,
                                          type = c("difference", "ratio"),
                                          nsim = 1,## large numbers need much more time
                                          sample_size = NULL,
                                          sample_frac = NULL,
                                          verbose = FALSE,
                                          progress = "none",
                                          parallel = TRUE,
                                          paropts = NULL)
    
    save(CNN_sim_impNULL,file = paste0("sim_",j,"_",para[i],"_model_CNN_L_Dropout_sim_impNULL.RData") )
    write.csv(CNN_sim_impNULL$CNN_Decrease_acc,file =paste0("sim_",j,"_",para[i],"_model_CNN_L_Dropout_impNULL_Decrease_acc.csv") )
    Sys.time()
    print(paste0("sim_",j,"_",para[i],"CNN optimal model permuting & NULL importance finished"))
    #registerDoSEQ()
    #stopCluster(cl)
  }
}


Sys.time()

save.image(file = "CNN_final_after_trained_data.RData")

print("ALl finished")
