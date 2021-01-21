library("caret")### for ML calling functions and performance estimation
library("caretEnsemble")
library("NeuralNetTools")
library("neuralnet")
library("DeepGenomeScan")

source("DLModels.R")
set.seed(123)
###pre-trial to figure out which type of model suits best for the data

### see DL_trials

### after selected a better model from model list, then check the optimal hyperparameters from it

### modify and change the search range (number of neurons changed to 0-20) according to the trials
### Note that this model used the developemnet version of neuralnet
k=2
a=0.8
linear=function(x) a*x
customRelu =function(x) {x/(1+exp(-2*k*x))} 
softplus =function(x) log(1 + exp(x)) 


mlpneuralnet1<- list(label = "Neural Network",
                      library = "neuralnet",
                      loop = NULL,
                      type = c('Regression'),
                      parameters = data.frame(parameter = c('layer1', 'layer2', 'layer3',"activation1","activation2","linear.output"),
                                              class = c('numeric', 'numeric', 'numeric',"character","character","character"),
                                              label = c('Number of hidden Units in Layer 1', 'number of hidden Units in Layer 2', 'number of hidden Units in Layer 3',"Activation function in hidden layer","Activation function in output layer","Activation function linear out choice")),
                      grid = function(x, y, len = NULL, search = "grid") {
                        afuncs=c("logistic", "tanh","softplus")
                        outputf=c("TRUE","FALSE")
                        if(search == "grid") {
                          out <- expand.grid(layer1 = ((1:len) * 2) - 1, layer2 = 0, layer3 = 0, activation1=c("logistic", "tanh","softplus"),activation2=c("logistic", "tanh","softplus"),linear.output=c("TRUE","FALSE"))
                        } else {
                          out <- data.frame(layer1 = sample(2:20, replace = TRUE, size = len),
                                            layer2 = sample(c(0, 2:20), replace = TRUE, size = len),
                                            layer3 = sample(c(0, 2:20), replace = TRUE, size = len),
                                            activation1=sample(afuncs, size = len, replace = TRUE),
                                            activation2=sample(afuncs, size = len, replace = TRUE),
                                            linear.output=sample(outputf,size = len,replace = TRUE))
                        }
                        out
                      },
                      fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                        colNames <- colnames(x)
                        dat <- if(is.data.frame(x)) x else as.data.frame(x, stringsAsFactors = TRUE)
                        dat$.outcome <- y
                        form <- as.formula(paste(".outcome ~",paste(colNames, collapse = "+")))
                        if(param$layer1 == 0) stop("the first layer must have at least one hidden unit")
                        if(param$layer2 == 0 & param$layer2 > 0) stop("the second layer must have at least one hidden unit if a third layer is specified")
                        nodes <- c(param$layer1)
                        if(param$layer2 > 0) {
                          nodes <- c(nodes, param$layer2)
                          if(param$layer3 > 0) nodes <- c(nodes, param$layer3)
                        }
                        actf=(param$activation1)### note the difference in c(param$activation) for this model and other model, becaue the self-defined softplus function can't be a vector, so here we should not use c(,softplus)
                        outputactf=(param$activation2)
                        linear.output=c(param$linear.output)
                        neuralnet::neuralnet(form, algorithm="rprop+",data = dat,rep=100, hidden = nodes, stepmax = 1e+09, learningrate.factor = list(minus = 0.5,plus = 1.2),act.fct=actf,output.act.fct=outputactf,linear.output=linear.output,...)
                      },## if the cluster do not use parallel and runs slowly, lower the number of repeat. 
                      predict = function(modelFit, newdata, submodels = NULL) {
                        newdata <- newdata[, modelFit$model.list$variables, drop = FALSE]
                        neuralnet::predict(modelFit, covariate = newdata)$net.result[,1]
                      },
                 varImp = function(object, ...){
                   imps <- NeuralNetTools::olden(object,bar_plot =FALSE)
                   out <- data.frame(Overall = as.vector(imps))
                   rownames(out) <- names(imps)
                   out
                 },
                      prob = NULL,
                      tags = c("Neural Network"),
                      sort = function(x) x[order(x$layer1, x$layer2, x$layer3,x$activation1,x$activation2,x$linear.output),])






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


