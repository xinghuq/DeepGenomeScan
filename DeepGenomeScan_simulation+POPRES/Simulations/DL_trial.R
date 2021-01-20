library("caret")
library("DeepGenomeScan")
library("caretEnsemble")
library("NeuralNetTools")
library("neuralnet")
library("h2o")
library("FCNN4R")
library("keras")
library("kerasR")
library("tensorflow")


source("DLModels.R")
set.seed(123)
###pre-trial to figure out which type of model suits best for the data

econtrol <- trainControl(## 5-fold CV, repeat 5 times
  method = "adaptive_cv",
  number = 5,
  adaptive = list(min = 5, alpha = 0.05,method = "gls", complete = TRUE),
  repeats = 5,search = "random")


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

#### pre-trial for seach the optimal model list for this data
## use a single simulated data to search simdat[[1]]
for(j in 1:length(modelist)) {
  for (i in 1:length(para)){
    print(paste0("model_",j,"_",para[i],"General search for model list"))
    Sys.time()
    #env=normalize(simdata[[j]][,2:11])
    genotype_norm=as.data.frame(apply(simdata[[1]],2,normalize))
    simf=as.formula(paste(colnames(env)[i],paste(names(genotype_norm[,15:1014]), collapse="+"),sep="~"))
    
    model_envi_mlp=DeepGenomeScan(simf, data=(genotype_norm),
                                method=modelist[[j]],
                                metric = "MAE",## "Accuracy", "RMSE"
                                preProcess=c("scale"),
                                tuneLength = 100, ### search 100 combinations of parameters
                                verbose=0,# verbose=1 is reporting the progress,o is sclience
                                trControl = econtrol)
    print(paste0("model_",j,"_",para[i],"tuning neuralnet finished"))
    Sys.time()
    save(model_envi_mlp,file=paste0("model_",j,"_",para[i],"_Scan_env_trained_model.RData"))
    ### olden importance
    neuralnet_simimp1=NeuralNetTools::olden(model_envi_mlp, bar_plot=FALSE)## This will ignore keras model
    write.csv(neuralnet_simimp1$scaled_importance,file = paste0("model_",j,"_",para[i],"_Scan_mlp_importance_env.csv"))
    write.csv(model_neuralnet_envi_mlp$results,file = paste0("model_",j,"_",para[i],"_model_mlp_envi_tuning.csv") )
    #registerDoSEQ()
    #stopCluster(cl)
  }
}

### compare model performance 
###"modelRSNNSmlpdecay"   "mlpFCNN4Rsgd"         "mlpneuralnet1"        "mlph2o"               "modelmlpkerasdropout"

load("model_modelRSNNSmlpdecay_env1_Scan_env_trained_model.RData")
Model_RSNNS=model_envi_mlp
load("model_mlpFCNN4Rsgd_env1_Scan_env_trained_model.RData")
Model_FCNN4=model_envi_mlp
load("model_mlpneuralnet1_env1_Scan_env_trained_model.RData")
Model_neuralnet=model_envi_mlp
load("model_mlph2o_env1_Scan_env_trained_model.RData")
Model_h2o=model_envi_mlp
load("model_modelmlpkerasdropout_env1_Scan_env_trained_model.RData")
Model_keras=model_envi_mlp


results <- resamples(list(Model_RSNNS=Model_RSNNS, Model_FCNN4=Model_FCNN4, Model_h2o=Model_h2o,Model_neuralnet=Model_neuralnet, Model_keras=Model_keras))

# boxplots of results
bwplot(results)
# dot plots of results
dotplot(results)
### do iteration do all the env factors to see which is best overall
