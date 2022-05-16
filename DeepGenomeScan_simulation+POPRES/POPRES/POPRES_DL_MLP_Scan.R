library(caret)
library(DeepGenomeScan)
library("caretEnsemble")
library("NeuralNetTools")
library("neuralnet")

set.seed(123)

#### perform deep genome scan on POPRES data

Sys.time()


normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

print("Loaing data")
Sys.time()

load("POPRES_DL_Scan_Data_KLFDAPC1.RData")
source("DLModels.R")

Sys.time()
print("setting DL_MLP")
########## MLP random search############################################
set.seed(123)
econtrol <- trainControl(## 5-fold CV, repeat 5 times
  method = "adaptive_cv",
  number = 5,
  ## repeated 5 times
  repeats = 5,
adaptive = list(min = 5, alpha = 0.05, method = "gls", complete = TRUE),
search = "random")

print("Start traing DL mlp")
Sys.time()
### loci related to long
print("DL_MLP regression on long")
Sys.time()
### I have changed the normalization process to the tarining process, 
#POPRES_DL_genotype_mat_norm=as.data.frame(POPRES_DL_genotype_mat$genotype)
genotype_dat=as.data.frame(POPRES_DL_genotype_mat$genotype)
colnames(genotype_dat)=POPRES_DL_genotype_mat$snp.id
rm(POPRES_DL_genotype_mat)
genotype_impute=caret::preProcess(genotype_dat,method = c("medianImpute")) ## impute the missing values
POPRES_DL_genotype_mat_norm=predict(genotype_impute, genotype_dat)

POPRES_mlpneuralnet_long<- DeepGenomeScan(x=POPRES_DL_genotype_mat_norm,y=popres_rm_I_ID$longitude,
                                       method=c("mlpFCNN4Rsgd"),
                                       metric = "MAE",## "Accuracy", "RMSE","MAE","R squred"
                                       tuneLength = 100, ### search 100 combinations of parameters
                                      # preProcess=c("knnImpute"), 
#verbose=0,# verbose=1 is reporting the progress,o is sclience
                                       trControl = econtrol)
    
    #  keras::unserialize_model(gbmGrid_model_mlpneuralnet_env$finalModel$object)
    POPRES_MLPImp_long <- as.data.frame(NeuralNetTools::olden(POPRES_mlpneuralnet_long, bar_plot=FALSE), stringsAsFactors = TRUE)
    write.csv(POPRES_MLPImp_long$importance,file = "POPRES_mlpneuralnet_long.csv")
    save(POPRES_mlpneuralnet_long,file="POPRES_mlpneuralnet_long.RData")
    POPRES_MLPImp_long2=OldenWeights_FCNN4R(POPRES_mlpneuralnet_long$finalModel)
    write.csv(POPRES_MLPImp_long2,file = paste0("sim_",j,"_",para[i],"_mlp_olden2_importance_long.csv"))
    write.csv(POPRES_mlpneuralnet_long$results,file = "POPRES_mlpneuralnet_long_tuning.csv")
     rm(POPRES_mlpneuralnet_long)
#### loci related to latitute
    print("DL_MLP regression on long finished")
    Sys.time()
    print("DL_MLP regression on lat ")
    Sys.time()
    
    POPRES_mlpneuralnet_lat<- DeepGenomeScan(x=POPRES_DL_genotype_mat_norm,y=popres_rm_I_ID$latitude,
                                        method=c("mlpFCNN4Rsgd"),
                                        metric = "MAE",## "Accuracy", "RMSE"
                                        tuneLength = 100, ### search > 100 combinations of parameters
                                       # preProcess=c("knnImpute"),
#verbose=0,# verbose=1 is reporting the progress,o is sclience
                                        trControl = econtrol)
  
    POPRES_MLPImp_lat <- as.data.frame(NeuralNetTools::olden(POPRES_mlpneuralnet_lat, bar_plot=FALSE), stringsAsFactors = TRUE)
    write.csv(POPRES_MLPImp_lat$importance,file = "POPRES_mlpneuralnet_lat.csv")
    save(POPRES_mlpneuralnet_lat,file="POPRES_mlpneuralnet_lat.RData")
    POPRES_MLPImp_lat2=OldenWeights_FCNN4R(POPRES_mlpneuralnet_lat$finalModel)
    write.csv(POPRES_MLPImp_lat2,file = paste0("sim_",j,"_",para[i],"_mlp_olden2_importance_lat.csv"))
    write.csv(POPRES_mlpneuralnet_lat$results,file = "POPRES_mlpneuralnet_lat_tuning.csv")
    print("DL_MLP regression on  lat finished")
    Sys.time() 
#stopCluster(cl)
rm(POPRES_mlpneuralnet_lat)
Sys.time()
print("DL_MLP regression on long and lat finished")
Sys.time()


save.image(file = "POPRES_MLP_Genome_Scan_true_location_trained_data.RData")


print("Loaing data")
Sys.time()

load("POPRES_DL_Scan_Data_KLFDAPC1.RData")


Sys.time()
print("setting DL_MLP")
########## MLP random search############################################
set.seed(123)
econtrol <- trainControl(## 5-fold CV, repeat 5 times
  method = "adaptive_cv",
  number = 5,
  ## repeated 5 times
  repeats = 5,
  adaptive = list(min = 5, alpha = 0.05, method = "gls", complete = TRUE),
  search = "random")

print("Start traing DL mlp")
Sys.time()
### loci related to long
print("DL_MLP regression on long")
Sys.time()
### i have changed the normalization process to the tarining process, 
#POPRES_DL_genotype_mat_norm=as.data.frame(POPRES_DL_genotype_mat$genotype)
genotype_dat=as.data.frame(POPRES_DL_genotype_mat$genotype)
colnames(genotype_dat)=POPRES_DL_genotype_mat$snp.id
rm(POPRES_DL_genotype_mat)
genotype_impute=caret::preProcess(genotype_dat,method = c("medianImpute"))
POPRES_DL_genotype_mat_norm=predict(genotype_impute, genotype_dat)

print("DL_MLP regression on KLFDAPC5 RD1-2")
print("DL_MLP regression on KLFDAPC5 RD1")
Sys.time()
POPRES_mlpneuralnet_KLFDAPC5_RD1<- DeepGenomeScan(x=POPRES_DL_genotype_mat_norm, y=KLFDAPC_popresnew5$data$PC1,
                                                   method=c("mlpFCNN4Rsgd"),
                                                   metric = "MAE",## "Accuracy", "RMSE"
                                                   tuneLength = 100, ### search 100 combinations of parameters
                                                   # preProcess=c("knnImpute"), 
                                                   #verbose=0,# verbose=1 is reporting the progress,o is sclience
                                                   trControl = econtrol)
Sys.time()

POPRES_MLPImp_KLFDAPC5_RD1 <- as.data.frame(NeuralNetTools::olden(POPRES_mlpneuralnet_KLFDAPC5_RD1, bar_plot=FALSE), stringsAsFactors = TRUE)
write.csv(POPRES_MLPImp_KLFDAPC5_RD1$importance,file = "IMP_POPRES_mlpneuralnet_KLFDAPC5_RD1.csv")
save(POPRES_mlpneuralnet_KLFDAPC5_RD1,file="POPRES_mlpneuralnet_KLFDAPC5_RD1.RData")
POPRES_MLPImp_KLFDAPC5_RD1_olden=OldenWeights_FCNN4R(POPRES_mlpneuralnet_KLFDAPC5_RD1$finalModel)
write.csv(POPRES_MLPImp_KLFDAPC5_RD1_olden,file = paste0("sim_",j,"_",para[i],"_IMP_POPRES_mlpneuralnet_KLFDAPC5_RD1_olden.csv"))

write.csv(POPRES_mlpneuralnet_KLFDAPC5_RD1$results,file = "POPRES_mlpneuralnet_KLFDAPC5_RD1_tuning.csv")

rm(POPRES_mlpneuralnet_KLFDAPC5_RD1)
print("DL_MLP regression on KLFDAPC5 RD2")
Sys.time()
#### RD2
POPRES_mlpneuralnet_KLFDAPC5_RD2<- DeepGenomeScan(x=POPRES_DL_genotype_mat_norm, y=KLFDAPC_popresnew5$data$PC2,
                                                   method=c("mlpFCNN4Rsgd"),
                                                   metric = "MAE",## "Accuracy", "RMSE"
                                                   tuneLength = 100, ### search 100 combinations of parameters
                                                   #   preProcess=c("knnImpute"),
                                                   #verbose=0,# verbose=1 is reporting the progress,o is sclience
                                                   trControl = econtrol)

#  keras::unserialize_model(gbmGrid_model_mlpneuralnet_env$finalModel$object)
POPRES_MLPImp_KLFDAPC5_RD2 <- as.data.frame(NeuralNetTools::olden(POPRES_mlpneuralnet_KLFDAPC5_RD2, bar_plot=FALSE), stringsAsFactors = TRUE)
write.csv(POPRES_MLPImp_KLFDAPC5_RD2$importance,file = "IMP_POPRES_mlpneuralnet_KLFDAPC5_RD2.csv")
save(POPRES_mlpneuralnet_KLFDAPC5_RD2,file="POPRES_mlpneuralnet_KLFDAPC5_RD2.RData")
POPRES_MLPImp_KLFDAPC5_RD2_olden=OldenWeights_FCNN4R(POPRES_mlpneuralnet_KLFDAPC5_RD2$finalModel)
write.csv(POPRES_MLPImp_KLFDAPC5_RD2_olden,file = paste0("sim_",j,"_",para[i],"_IMP_POPRES_mlpneuralnet_KLFDAPC5_RD2_olden.csv"))
write.csv(POPRES_mlpneuralnet_KLFDAPC5_RD2$results,file = "POPRES_mlpneuralnet_KLFDAPC5_RD2_tuning.csv")
Sys.time()
rm(POPRES_mlpneuralnet_KLFDAPC5_RD2)
print("DL_MLP regression on KLFDAPC5 RD1-2 finished")
Sys.time()
print("DL_MLP Genome scan all finished")

Sys.time()
save.image(file = "POPRES_MLP_Genome_Scan_all_true_location_KLFDAPC5_trained_data.RData")
Sys.time()

