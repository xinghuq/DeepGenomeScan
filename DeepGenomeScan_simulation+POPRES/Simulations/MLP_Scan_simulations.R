library("caret")### for ML calling functions and performance estimation
library("caretEnsemble")
library("NeuralNetTools")
library("neuralnet")
library("DeepGenomeScan")

source("DLModels.R")
set.seed(123)
###pre-trial to figure out which type of model suits best for the data

options(warn=-1)
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

dirs=list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
fls <- list.files(dirs, pattern="*.csv", full.names = TRUE)
simdata=lapply(fls, read.csv)
#save(simdata,file = "Dimulated_Data.RData")
#### separate data, while this is not necrssary as long as you know and subset the data
env=apply((simdata[[1]][,2:11]),2,normalize)
genotype=simdata[[1]][,-(1:14)]
qtrait=simdata[[1]][,12:14]
para=colnames(env)

### model pre-trial in DL_trial to select which model could performe better


#### Please note that we used development version of neuralnet https://github.com/bips-hb/neuralnet, which can specify the out.act.fct. 
## if you install neuralnet from CRAN, the source model neuralnet1 will not allow you to change output.act.fct

########## neuralnet mlp random search############################################
set.seed(123)
#econtrol <- caret::trainControl(## 5-fold CV, repeat 5 times
 # method = "adaptive_cv",
 # number = 5,
  ## repeated ten times
 # repeats = 5,
 # adaptive = list(min = 5, alpha = 0.05,method = "gls", complete = TRUE),
 # search = "random")


### we use repeatedcv for simulations 

econtrol <- trainControl(## 5-fold CV, repeat 5 times
  method = "repeatedcv",
  number = 5,
  repeats = 5,search = "random")

feature_names=colnames(genotype)
print("Start traing DL_MLP")
Sys.time()
softplus <- function(x) log(1 + exp(x))

for(j in 1:length(simdata)) {
  for (i in 1:length(para)){
    print(paste0("sim_",j,"_",para[i],"tuning RNN 100 parameters"))
    Sys.time()
    #env=normalize(simdata[[j]][,2:11])
    genotype_norm=as.data.frame(apply(simdata[[j]],2,normalize))
    simf=as.formula(paste(colnames(env)[i],paste(names(genotype_norm[,15:1014]), collapse="+"),sep="~"))
    
    model_neuralnet_envi_mlp=caret::train(simf, data=(genotype_norm),
                                    method=mlpneuralnet1,
                                    metric = "MAE",## "Accuracy", "RMSE","MAE","R-squared"
                                  #  preProcess=c("scale"),
                                    tuneLength = 100, ### search 100 combinations of parameters
                                    verbose=0,# verbose=1 is reporting the progress,o is sclience
                                    trControl = econtrol)
    print(paste0("sim_",j,"_",para[i],"tuning neuralnet finished"))
    Sys.time()
    save(model_neuralnet_envi_mlp,file=paste0("sim_",j,"_",para[i],"_neuralnet_Scan_env_trained_model.RData"))
    ### olden importance
    neuralnet_simimp1=NeuralNetTools::olden(model_neuralnet_envi_mlp, bar_plot=FALSE)
    write.csv(neuralnet_simimp1$importance,file = paste0("sim_",j,"_",para[i],"_mlpneuralnet_importance_env.csv"))
    write.csv(model_neuralnet_envi_mlp$results,file = paste0("sim_",j,"_",para[i],"_model_neuralnet_mlp_envi_tuning.csv") )
    #registerDoSEQ()
    #stopCluster(cl)
  }
}


Sys.time()

save.image(file = "neuralnet_mlp_Scan_sims_final_after_trained_data.RData")


