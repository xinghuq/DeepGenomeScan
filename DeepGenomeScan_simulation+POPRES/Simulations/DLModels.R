#### This source contains all alternative models constructed with different libraries for deep learning based genome scan in our DeepGenomeScan
modelist=c("modelRSNNSmlpdecay","mlpFCNN4Rsgd","mlpneuralnet1","mlph2o","modelmlpkerasdropout")
set.seed(123)
####modelRSNNSmlpdecay
 library(RSNNS)
modelRSNNSmlpdecay<- list(label = "Multi-Layer Perceptron, multiple layers",
                           library = "RSNNS",
                           loop = NULL,
                           type = c('Regression', 'Classification'),
                           parameters = data.frame(parameter = c('layer1','layer2','layer3', 'decay',"activation1","activation2"),
                                                   class = c('numeric','numeric','numeric', 'numeric',"character","character"),
                                                   label = c('number of hidden Units layer1','number of hidden Units layer2','number of hidden Units layer3', 'Weight Decay',"hiddenActFunc","outputActFunc")),
                           
                           grid = function(x, y, len = NULL, search = "grid"){
                             afuncs <- c("Act_Logistic","Act_TanH","Act_RBF_Gaussian","Act_IdentityPlusBias")
                             if(search == "grid") {
                               
                               out <- expand.grid(layer1 = ((1:len) * 2) - 1, layer2 = ((1:len) * 2) - 1, layer3 = ((1:len) * 2) - 1,
                                                  decay = c(0, 10 ^ seq(-1, -4, length = len - 1)),activation1=c("Act_Logistic","Act_TanH","Act_RBF_Gaussian","Act_IdentityPlusBias"),activation2=c("Act_Logistic","Act_TanH","Act_RBF_Gaussian","Act_IdentityPlusBias"))
                             } else {
                               out <- data.frame(layer1 = sample(2:100, replace = TRUE, size = len),
                                                 layer2 = sample(c(0, 2:100), replace = TRUE, size = len),
                                                 layer3 = sample(c(0, 2:100), replace = TRUE, size = len),
                                                 decay = 10^runif(len, min = -5, max = 1),
                                                 activation1 = sample(afuncs, size = len, replace = TRUE),
                                                 activation2 = sample(afuncs, size = len, replace = TRUE)
                               )
                             }
                             out
                           },
                           fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                             theDots <- list(...)
                             theDots <- theDots[!(names(theDots) %in% c("size","linOut"))]
                             if(any(names(theDots) == "learnFunc"))
                             {
                               theDots$learnFunc <- NULL
                               warning("Cannot over-ride 'learnFunc' argument for this model. BackpropWeightDecay is used.")
                             }
                             if(any(names(theDots) == "learnFuncParams"))
                             {
                               prms <- theDots$learnFuncParams
                               prms[2] <-  param$decay
                               warning("Over-riding weight decay value in the 'learnFuncParams' argument you passed in. Other values are retained")
                             } else prms <- c(0.2, param$decay, 0.0, 0.0)
                             
                             if(is.factor(y)) {
                               y <- RSNNS:::decodeClassLabels(y)
                               lin <- FALSE
                             } else lin <- TRUE
                             
                             nodes <- c(param$layer1, param$layer2, param$layer3)
                             if (any(nodes == 0)) {
                               nodes <- nodes[nodes > 0]
                               warning(
                                 "At least one layer had zero units and ",
                                 "were removed. The new structure is ",
                                 paste0(nodes, collapse = "->"), call. = FALSE
                               ) 
                             }
                             func1= c(param$activation1)
                             func2= c(param$activation2)
                             
                             args <- list(x = x,
                                          y = y,
                                          learnFunc = "BackpropWeightDecay",
                                          learnFuncParams = prms,
                                          hiddenActFunc = func1,
                                          size = nodes,
                                          outputActFunc = func2,
                                          linOut = lin)
                             args <- c(args, theDots)
                             do.call(RSNNS::mlp, args)
                           },
                           predict = function(modelFit, newdata, submodels = NULL) {
                             out <- predict(modelFit, newdata)
                             if(modelFit$problemType == "Classification")
                             {
                               out <- modelFit$obsLevels[apply(out, 1, which.max)]
                             } else out <- out[,1]
                             out
                           },
                  varImp = function(object, ...){
                    imps <- NeuralNetTools::olden(object,bar_plot =FALSE)
                    out <- data.frame(Overall = as.vector(imps))
                    rownames(out) <- names(imps)
                    out
                  },
                           prob = function(modelFit, newdata, submodels = NULL) {
                             out <- predict(modelFit, newdata)
                             colnames(out) <- modelFit$obsLevels
                             out
                           },
                           levels = function(x) x$obsLevels,
                           tags = c("Neural Network","L2 Regularization"),
                           sort = function(x) x[order(x$layer1, x$layer2, x$layer3, -x$decay),])

 
library(FCNN4R) #### should install from source, or alternative from devtools::install_github("xinghuq/Models/FCNN4R")
##mlpFCNN4Rsgd
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
                         out <- data.frame(units1 = sample(2:100, replace = TRUE, size = len),
                                           units2 = sample(0:100, replace = TRUE, size = len),
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

#### mlpneuralnet1
library(neuralnet)
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
                          out <- data.frame(layer1 = sample(2:100, replace = TRUE, size = len),
                                            layer2 = sample(c(0, 2:100), replace = TRUE, size = len),
                                            layer3 = sample(c(0, 2:100), replace = TRUE, size = len),
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
                        neuralnet::neuralnet(form, algorithm="rprop+",data = dat,rep=1, hidden = nodes, stepmax = 1e+09, learningrate.factor = list(minus = 0.5,plus = 1.2),act.fct=actf,output.act.fct=outputactf,linear.output=linear.output,...)
                      },
                      predict = function(modelFit, newdata, submodels = NULL) {
                        newdata <- newdata[, modelFit$model.list$variables, drop = FALSE]
                        predict(modelFit, covariate = newdata)$net.result[,1] ### neuralnet::predict,or neuralnet::compute for old version 
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


##mlph2o
mlph2o<- list(label = "DL_h2o",
                         library = "h2o",
                         type = c("Regression", "Classification"),
                         parameters = data.frame(parameter = c('units1','units2', 'l2reg', "rho", "activation"),
                                                 class = c(rep('numeric', 4), rep("character",1)),
                                                 label = c('Number of hidden Units1','Number of hidden Units2', 'L2 Regularization', 
                                                           "Adaptive learning rate time decay factor", 
                                                           'Activation function at hidden layer')),
                         grid = function(x, y, len = NULL, search = "grid") {
                           
                           afuncs=c("Tanh", "TanhWithDropout", "Rectifier", "RectifierWithDropout")
                           
                           if(search == "grid") {
                             out <- expand.grid(units1 = ((1:len) * 2) - 1, ### two hidden layers
                                                units2 = ((1:len) * 2) - 1, 
                                                l2reg = c(0, 10 ^ seq(-1, -4, length = len - 1)), 
                                                rho = 2e-3, 
                                                activation=c("Tanh", "TanhWithDropout", "Rectifier", "RectifierWithDropout"))
                           } else {
                             out <- data.frame(units1 = sample(1:100, replace = TRUE, size = len),
                                               units2 = sample(0:100, replace = TRUE, size = len),
                                               l2reg = 10^runif(len, min = -5, 1),
                                               rho = runif(len),
                                               activation= sample(afuncs, size = len, replace = TRUE))
                           }
                           out
                         },
                         loop = NULL,
                         fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                           
                           dat <- if(!is.data.frame(x)) as.data.frame(x, stringsAsFactors = TRUE) else x
                           dat$.outcome <- y
                           p <- ncol(dat)
                           frame_name <- paste0("tmp_DL_h2o_dat_",sample.int(100000, 1))
                           tmp_train_dat = h2o::as.h2o(dat, destination_frame = frame_name)
                           
                           nodes <- c(param$units1)
                           if(param$units2 > 0) {
                             nodes <- c(nodes, param$units2)
                           }
                           acts= as.character(param$activation)
                           out <- h2o::h2o.deeplearning(x = colnames(x), y = ".outcome",
                                                        training_frame = tmp_train_dat,
                                                        standardize = F,
                                                        # model_id = "deep_model",        
                                                        activation =acts, 
                                                        hidden=nodes,
                                                        epochs = 100,       
                                                        seed = 123,
                                                        variable_importances = T,
                                                        l2=param$L2reg,
                                                        stopping_metric=ifelse(is.factor(env$envir1), "misclassification", "RMSE"),
                                                        rho=param$rho,...)
                           h2o::h2o.getModel(out@model_id)
                         },
                         predict = function(modelFit, newdata, submodels = NULL) {
                           frame_name <- paste0("new_DL_h2o_dat_",sample.int(100000, 1))
                           newdata <- h2o::as.h2o(newdata, destination_frame = frame_name)
                           as.data.frame(h2o::h2o.predict(modelFit, newdata), stringsAsFactors = TRUE)[,1]
                         },
                         prob = function(modelFit, newdata, submodels = NULL) {
                           frame_name <- paste0("new_DL_h2o_dat_",sample.int(100000, 1))
                           newdata <- h2o::as.h2o(newdata, destination_frame = frame_name)
                           as.data.frame(h2o::h2o.predict(modelFit, newdata), stringsAsFactors = TRUE)[,-1]
                         },
                         predictors = function(object, ...) {
                           out <- as.data.frame(h2o::h2o.varimp(object$finalModel), stringsAsFactors = TRUE)
                           colnames(out)[colnames(out) == "scaled_importance"] <- "Overall"
                           out <- out[!is.na(out$Overall),]   
                           out$names
                         },
                         varImp = function(object, ...) {
                           out <- as.data.frame(h2o::h2o.varimp(object$finalModel), stringsAsFactors = TRUE)
                           colnames(out)[colnames(out) == "scaled_importance"] <- "Overall"
                           rownames(out) <- out$names
                          # out <- out[!is.na(out$Overall), c("Overall"), drop = FALSE]   
                         #  all_var <- object$finalModel@allparameters$x
                          # if(any(!(all_var %in% rownames(out)))) {
                         #    missing <- all_var[!(all_var %in% rownames(out))]
                        #    tmp <- data.frame(OVerall = rep(0, length(missing)))
                         #    rownames(tmp) <- missing
                         #    out <- rbind(out, tmp)
                          # }
                           out
                         },
                         levels = NULL,
                         tags = c("Deep Learning", "MLP", 
                                  "L2 Regularization", "learning rate",
                                  "neural network regression"),
                         sort = function(x) x[order(-x$units1, x$units2),],
                         trim = NULL)


#modelmlpkerasdropout
library(keras)
library(kerasR)
modelmlpkerasdropout<- list(label = "Multilayer Perceptron Network",
                             library = "keras",
                             loop = NULL,
                             type = c('Regression', "Classification"),
                             parameters = data.frame(
                               parameter = c('units1',"units2", 'dropout1',"dropout2","batch_size",
                                             "lr", "rho", "decay", "activation1","activation2",
                                             "activation3"),
                               class = c(rep('numeric', 8), rep("character",3)),
                               label = c('Number of hidden Units1', 'Number of hidden Units2','Dropout Rate1','Dropout Rate2', 
                                         "Batch Size", "Learning Rate",
                                         "Rho", "Learning Rate Decay",
                                         "Activation Function1","Activation Function2","Activation Function3")
                             ),
                             grid = function(x, y, len = NULL, search = "grid") {
                               afuncs <- c("sigmoid", "relu", "tanh")
                               if(search == "grid") {
                                 out <- expand.grid(
                                   units1 = ((1:len) * 5) - 1, 
                                   units2=((1:len) * 5) - 1, 
                                   dropout1 = seq(0, .5, length = len), ### dropout rate will give warning if it is larger than 0.5
                                   dropout1 = seq(0, .5, length = len),
                                   batch_size = floor(nrow(x)/3),
                                   lr = 2e-6,
                                   rho = .8,
                                   decay = 0,
                                   activation1 = c("sigmoid", "relu", "tanh"), ### softmax is usually used for classification
                                   activation2 = c("sigmoid", "relu", "tanh"),
                                   activation3 = c("sigmoid", "relu", "tanh")
                                 )
                               } else {
                                 n <- nrow(x)
                                 out <- data.frame(
                                   units1 = sample(2:100, replace = TRUE, size = len), ### can be 0:100, so that can be 0 in a layer
                                   units2 = sample(0:100, replace = TRUE, size = len),### can be 0:100, so that can be 0 in a layer
                                   ## if you need, set the last layer as units3: units3 = sample(2:25, replace = TRUE, size = len),
                                   dropout1 = runif(len, max = .5),
                                   dropout2 = runif(len, max = .5), 
                                   batch_size = floor(n*runif(len, min = .1)),
                                   lr = runif(len),
                                   rho = runif(len),
                                   decay = 10^runif(len, min = -5, 0),
                                   activation1 = sample(afuncs, size = len, replace = TRUE),
                                   activation2 = sample(afuncs, size = len, replace = TRUE),
                                   activation3 = sample(afuncs, size = len, replace = TRUE)
                                 )
                               }
                               out
                             },
                             
                             ### construct MLP model
                             
                             fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                               require(dplyr)
                               K <- keras::backend()
                               K$clear_session()
                               if(!is.matrix(x)) x <- as.matrix(x)
                               model=keras_model_sequential() 
                               model %>%
                                 layer_dense(units = param$units1, activation = as.character(param$activation1), input_shape =ncol(x)) %>% 
                                 layer_dropout(rate = param$dropout1,seed = sample.int(1000, 1)) %>% 
                                 layer_dense(units = param$units2, activation = as.character(param$activation2))%>%
                                 layer_dropout(rate = param$dropout2,seed = sample.int(1000, 1)) 
                               
                               if(is.factor(y)) {
                                 y <- class2ind(y)
                                 model %>% 
                                   keras::layer_dense(
                                     units = length(lev), 
                                     activation = 'softmax'
                                   ) %>%
                                   keras::compile(
                                     loss = "categorical_crossentropy",
                                     optimizer = keras::optimizer_rmsprop(
                                       lr = param$lr,
                                       rho = param$rho,
                                       decay = param$decay
                                     ),
                                     metrics = "accuracy"
                                   )
                               } else {
                                 model %>% 
                                   keras::layer_dense(units = 1, activation = as.character(param$activation3)) %>%
                                   keras::compile(
                                     loss = "mean_squared_error",
                                     optimizer = keras::optimizer_rmsprop(
                                       lr = param$lr,
                                       rho = param$rho,
                                       decay = param$decay
                                     ),
                                     metrics = "mean_squared_error"
                                   )
                               }
                               ## alternative using sgd  #Compile the model
                               #     model%>% compile(
                               #    optimizer='sgd',
                               #    loss='mean_squared_error',
                               #     metrics='mean_squared_error')
                               
                               
                               model %>% keras::fit(
                                 x = x, 
                                 y = y,
                                 batch_size = param$batch_size,
                                 ...
                               )
                               if(last)
                                 model <- keras::serialize_model(model)
                               list(object = model)
                             },
                             predict = function(modelFit, newdata, submodels = NULL) {
                               if(inherits(modelFit$object, "raw"))
                                 modelFit$object <- keras::unserialize_model(modelFit$object)
                               if(!is.matrix(newdata)) 
                                 newdata <- as.matrix(newdata)
                               out <- predict(modelFit$object, newdata)
                               ## check for model type
                               if(ncol(out) == 1) {
                                 out <- out[, 1]
                               } else {
                                 out <- modelFit$obsLevels[apply(out, 1, which.max)]
                               }
                               out
                             },
                             prob =  function(modelFit, newdata, submodels = NULL) {
                               if(inherits(modelFit$object, "raw"))
                                 modelFit$object <- keras::unserialize_model(modelFit$object)
                               if(!is.matrix(newdata)) 
                                 newdata <- as.matrix(newdata)
                               out <- predict(modelFit$object, newdata)
                               colnames(out) <- modelFit$obsLevels
                               as.data.frame(out, stringsAsFactors = TRUE)
                             },
                  
                  varIMP=NULL,
  #                feature_names=colnames(genotype),     ### see example, we could not expect the var-importance hang the time here, but it's better used the final model to estimete 
 #                 varImp = function(object, ...) {
   #                 out <- as.data.frame(MLP_varIMP_NULL_model(object,
  #                                                             feature_names,
  #                                                             train_y=env,
   #                                                            train_x=as.matrix(genotype),
 #                                                              #smaller_is_better = FALSE,
  #                                                             type = c("difference", "ratio"),
 #                                                              nsim = 100,## MCMC permutation large numbers need much more time
 #                                                              sample_size = NULL,
 #                                                              sample_frac = NULL,
  #                                                             verbose = FALSE,
  ##                                                             progress = "none",
   #                                                            parallel = TRUE,
  #                                                             paropts = NULL))
  #                  colnames(out$MLP_Decrease_acc)[colnames(out$MLP_Decrease_acc) == "relative_importance"] <- "Overall"
 #                   rownames(out$MLP_Decrease_acc) <- out$variable
 #                   out[, c("Overall"), drop = FALSE]
#                   }
                  
                             tags = c("Neural Network"),
                             sort = function(x) x[order(x$units1, x$units2,x$dropout1,x$dropout2),],
                             notes = paste("After `train` completes, the keras model object is serialized",
                                           "so that it can be used between R session. When predicting, the", 
                                           "code will temporarily unsearalize the object. To make the", 
                                           "predictions more efficient, the user might want to use ", 
                                           "`keras::unsearlize_model(object$finalModel$object)` in the current", 
                                           "R session so that that operation is only done once.",
                                           "Also, this model cannot be run in parallel due to",
                                           "the nature of how tensorflow does the computations.",
                                           
                                           "Unlike other packages used by `train`, the `dplyr`",
                                           "package is fully loaded when this model is used."),
                             check = function(pkg) {
                               testmod <- try(keras::keras_model_sequential(),
                                              silent = TRUE)
                               if(inherits(testmod, "try-error"))
                                 stop("Could not start a sequential model. ",
                                      "`tensorflow` might not be installed. ",
                                      "See `?install_tensorflow`.", 
                                      call. = FALSE)
                               TRUE
                             })
