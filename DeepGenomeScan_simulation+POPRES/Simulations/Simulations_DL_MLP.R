library("caret")### for ML calling functions and performance estimation
library("caretEnsemble")
library("NeuralNetTools")
# library("FCNN4R") has to be installed
library("DeepGenomeScan")

source("DLModels.R")
set.seed(123)
###pre-trial to figure out which type of model suits best for the data

### see DL_trials

### After selecting a better model from model list, then check the optimal hyperparameters from it

### Modifying and changing the search range according to the trials
### 

mlpFCNN4Rsgd <- list(label = "Multilayer Perceptron Network by Stochastic Gradient Descent",
                     library = c("FCNN4R", "plyr"),
                     loop = NULL,
                     type = c('Regression', "Classification"),
                     parameters = data.frame(parameter = c('units1','units2', 'l2reg', 'lambda', "learn_rate", 
                                                           "momentum", "gamma", "minibatchsz", "repeats","activation1","activation2"),
                                             class = c(rep('numeric', 9), rep("character",2)),
                                             label = c('Number of hidden Units1','Number of hidden Units2', 'L2 Regularization', 
                                                       'RMSE Gradient Scaling', "Learning Rate", 
                                                       "Momentum", "Learning Rate Decay", "Batch Size",
                                                       "#Models",'Activation function at hidden layer','Activation function at output layer')),
                     grid = function(x, y, len = NULL, search = "grid") {
                       n <- nrow(x)
                       afuncs=c("linear", "sigmoid", "tanh")
                       if(search == "grid") {
                         out <- expand.grid(units1 = ((1:len) * 2) - 1, 
                                            units2 = ((1:len) * 2) - 1, 
                                            l2reg = c(0, 10 ^ seq(-1, -4, length = len - 1)), 
                                            lambda = 0,
                                            learn_rate = 2e-6, 
                                            momentum = 0.9, 
                                            gamma = seq(0, .9, length = len),
                                            minibatchsz = floor(nrow(x)/3),
                                            repeats = 1,
                                            activation1=c("linear", "sigmoid", "tanh"),
                                            activation2=c("linear", "sigmoid", "tanh"))
                       } else {
                         out <- data.frame(units1 = sample(2:20, replace = TRUE, size = len),
                                           units2 = sample(0:20, replace = TRUE, size = len),
                                           l2reg = 10^runif(len, min = -5, 1),
                                           lambda = runif(len, max = .4),
                                           learn_rate = runif(len),
                                           momentum = runif(len, min = .5),
                                           gamma = runif(len),
                                           minibatchsz = floor(n*runif(len, min = .1)),
                                           repeats = sample(1:10, replace = TRUE, size = len),
                                           activation1 = sample(afuncs, size = len, replace = TRUE),
                                           activation2 = sample(afuncs, size = len, replace = TRUE))
                       }
                       out
                     },
                     fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                       if(!is.matrix(x)) x <- as.matrix(x)
                       if(is.factor(y)) {
                         y <- class2ind(y)
                         net <- FCNN4R::mlp_net(c(ncol(x), param$units1, param$units2,ncol(y)))
                         net <- FCNN4R::mlp_set_activation(net, layer = "h", activation = param$activation1)
                         net <- FCNN4R::mlp_set_activation(net, layer = "o", activation = param$activation2)
                         
                       } else {
                         y <- matrix(y, ncol = 1)
                         net <- FCNN4R::mlp_net(c(ncol(x), param$units1,param$units2, 1))
                         net <- FCNN4R::mlp_set_activation(net, layer = "h", activation = param$activation1)
                         net <- FCNN4R::mlp_set_activation(net, layer = "o", activation =param$activation2)
                       }
                       n <- nrow(x)
                       args <- list(net = net, 
                                    input = x, output = y, 
                                    learn_rate = param$learn_rate,
                                    minibatchsz = param$minibatchsz,
                                    l2reg = param$l2reg,
                                    lambda = param$lambda,
                                    gamma = param$gamma,
                                    momentum = param$momentum)
                       the_dots <- list(...) 
                       if(!any(names(the_dots) == "tol_level")) {
                         if(ncol(y) == 1) 
                           args$tol_level <- sd(y[,1])/sqrt(nrow(y)) else
                             args$tol_level <- .001
                       } 
                       
                       if(!any(names(the_dots) == "max_epochs")) 
                         args$max_epochs <- 1000
                       args <- c(args, the_dots)
                       out <- list(models = vector(mode = "list", length = param$repeats))
                       for(i in 1:param$repeats) {
                         args$net <- FCNN4R::mlp_rnd_weights(args$net)
                         out$models[[i]] <- do.call(FCNN4R::mlp_teach_sgd, args)
                       }
                       out
                     },
                     predict = function(modelFit, newdata, submodels = NULL) {
                       if(!is.matrix(newdata)) newdata <- as.matrix(newdata)
                       out <- lapply(modelFit$models, 
                                     function(obj, newdata)
                                       FCNN4R::mlp_eval(obj$net, input = newdata),
                                     newdata = newdata)
                       if(modelFit$problemType == "Classification") {
                         out <- as.data.frame(do.call("rbind", out), stringsAsFactors = TRUE)
                         out$sample <- rep(1:nrow(newdata), length(modelFit$models))
                         out <- plyr::ddply(out, plyr::`.`(sample), function(x) colMeans(x[, -ncol(x)]))[, -1]
                         out <- modelFit$obsLevels[apply(out, 1, which.max)]
                       } else {
                         out <- if(length(out) == 1) 
                           out[[1]][,1]  else {
                             out <- do.call("cbind", out)
                             out <- apply(out, 1, mean)
                           }
                       }
                       out
                     },
                     prob =  function(modelFit, newdata, submodels = NULL) {
                       if(!is.matrix(newdata)) newdata <- as.matrix(newdata)
                       out <- lapply(modelFit$models, 
                                     function(obj, newdata)
                                       FCNN4R::mlp_eval(obj$net, input = newdata),
                                     newdata = newdata)
                       out <- as.data.frame(do.call("rbind", out), stringsAsFactors = TRUE)
                       out$sample <- rep(1:nrow(newdata), length(modelFit$models))
                       out <- plyr::ddply(out, plyr::`.`(sample), function(x) colMeans(x[, -ncol(x)]))[, -1]
                       out <- t(apply(out, 1, function(x) exp(x)/sum(exp(x))))
                       colnames(out) <- modelFit$obsLevels
                       as.data.frame(out, stringsAsFactors = TRUE)
                     },
                     varImp = function(object, ...) {
                       imps <- lapply(object$models, caret:::GarsonWeights_FCNN4R, xnames = object$xNames)
                       imps <- do.call("rbind", imps)
                       imps <- apply(imps, 1, mean, na.rm = TRUE)
                       imps <- data.frame(var = names(imps), imp = imps)
                       imps <- plyr::ddply(imps, plyr::`.`(var), function(x) c(Overall = mean(x$imp)))
                       rownames(imps) <- as.character(imps$var)
                       imps$var <- NULL
                       imps[object$xNames,,drop = FALSE]
                     },
                     tags = c("Neural Network", "L2 Regularization"),
                     sort = function(x) x[order(x$units1, x$units2,-x$l2reg, -x$gamma),])

varImpFCNN4 = function(object, ...) {
  imps <- lapply(object$models, caret:::GarsonWeights_FCNN4R, xnames = object$xNames)
  imps <- do.call("rbind", imps)
  imps <- apply(imps, 1, mean, na.rm = TRUE)
  imps <- data.frame(var = names(imps), imp = imps)
  imps <- plyr::ddply(imps, plyr::`.`(var), function(x) c(Overall = mean(x$imp)))
  rownames(imps) <- as.character(imps$var)
  imps$var <- NULL
  imps[object$xNames,,drop = FALSE]
}

OldenWeights_FCNN4R=function (object, xnames = NULL, ynames = NULL) 
{
  beta <- (object$net@m_w_values[which(object$net@m_w_flags != 0L)])
  dims <- object$net@m_layers
  index <- (dims[1] + 1) * dims[2]
  i2h <- t(matrix(beta[1:index], ncol = dims[2]))
  i2h <- i2h[, -1, drop = FALSE]
  h2o <- matrix(beta[(index + 1):length(beta)], ncol = dims[3])
  h2o <- h2o[-1, , drop = FALSE]
  imp <- matrix(NA, nrow = dims[1], ncol = dims[3])
  for (output in 1:dims[3]) {
    Pij <- i2h * NA
    for (hidden in 1:dims[2]) Pij[hidden, ] <- i2h[hidden, ] * h2o[hidden, output] ### connected weights
#    Qij <- Pij * NA
 #   for (hidden in 1:dims[2]) Qij[hidden, ] <- Pij[hidden, ]/sum(Pij[hidden, ])

    Olden_sum=apply(Pij, 2, sum)
   # Sj <- apply(Qij, 2, sum)##Garson  imp[, output] <- Sj/sum(Sj) * 100
    imp[, output] <- Olden_sum/sum(Olden_sum) * 100 ## the relative importance X/sum(X)*100
#    rm(Pij, Qij, Sj)
  }
  rownames(imp) <- if (is.null(xnames)) 
    paste("X", 1:dims[1], sep = "")
  else xnames
  colnames(imp) <- if (is.null(ynames)) 
    paste("Y", 1:dims[3], sep = "")
  else ynames
  imp
}


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
#genotype=simdata[[1]][,-(1:14)]
#qtrait=simdata[[1]][,12:14]
para=colnames(env)

### model pre-trial in DL_trial to select which model could performe better

########## MLP random search############################################
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
    env=normalize(simdata[[j]][,2:11])
    genotype_norm=as.data.frame(apply(simdata[[j]],2,normalize))
    # simf=as.formula(paste(colnames(env)[i],paste(names(genotype_norm[,15:1014]), collapse="+"),sep="~"))
    
    model_neuralnet_envi_mlp=caret::train(y=env[,para[i]], x=(genotype_norm),
                                    method="mlpFCNN4Rsgd", ## or using "mlpSGD","mlph2o", the best practice are these FCNN4R and h2o
                                    metric = "RMSE",## "Accuracy", "RMSE","MAE","R-squared"
                                  #  preProcess=c("scale"),
                                    tuneLength = 100, ### search 100 combinations of parameters
                                    #verbose=0,# verbose=1 is reporting the progress,o is sclience
                                    trControl = econtrol)
    print(paste0("sim_",j,"_",para[i],"tuning MLP finished"))
    Sys.time()
    save(model_neuralnet_envi_mlp,file=paste0("sim_",j,"_",para[i],"_MLP_Scan_env_trained_model.RData"))
    write.csv(model_neuralnet_envi_mlp$results,file = paste0("sim_",j,"_",para[i],"_model_MLP_mlp_envi_tuning.csv") )  
    
    garson_simimp=varImpFCNN4(model_neuralnet_envi_mlp, bar_plot=FALSE)
    write.csv(garson_simimp$importance,file = paste0("sim_",j,"_",para[i],"_mlp_garson_importance_env.csv"))
    ### olden importance
    neuralnet_simimp1=NeuralNetTools::olden(model_neuralnet_envi_mlp, bar_plot=FALSE)
    write.csv(neuralnet_simimp1$importance,file = paste0("sim_",j,"_",para[i],"_mlp_olden1_importance_env.csv"))
    neuralnet_simimp2=OldenWeights_FCNN4R(model_neuralnet_envi_mlp$finalModel)
    write.csv(neuralnet_simimp2$importance,file = paste0("sim_",j,"_",para[i],"_mlp_olden2_importance_env.csv"))
 
    

    #registerDoSEQ()
    #stopCluster(cl)
  }
}


Sys.time()

save.image(file = "neuralnet_mlp_Scan_sims_final_after_trained_data.RData")


