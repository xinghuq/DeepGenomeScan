########################### compile with the framework
#modelmlpkerasdropout

 modelInfo<- list(label = "Multilayer Perceptron Network",
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
                                   units1 = sample(2:20, replace = TRUE, size = len), ### can be 0:100, so that can be 0 in a layer
                                   units2 = sample(2:20, replace = TRUE, size = len),### can be 0:100, so that can be 0 in a layer
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