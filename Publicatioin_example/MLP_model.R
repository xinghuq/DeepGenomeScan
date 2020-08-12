#### define the MLP model
nSNP=ncol(simdata[[1]][,-(1:14)]) 

modelmlpkerasdropout <- list(label = "Multilayer Perceptron Network",
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
                        units1 = sample(2:100, replace = TRUE, size = len),
                        units2 = sample(2:50, replace = TRUE, size = len),
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
                    model <- keras::keras_model_sequential()
                    model %>%   
                      keras_model_sequential() %>%
                    layer_dense(units = param$units1, activation = param$activation1, input_shape = c(nSNP)) %>% 
                    layer_dropout(rate = Dropout1,seed = sample.int(1000, 1)) %>% 
                    layer_dense(units = param$units2, activation = param$activation2) %>%
                    layer_dropout(rate = Dropout2,seed = sample.int(1000, 1)) 
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
                        keras::layer_dense(
                            layer_dense(units = 1, activation = param$activation3)
                        ) %>%
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
                  varImp = NULL,### see RNNSN
                  tags = c("Neural Network"),
                  sort = function(x) x[order(x$size, -x$dropout),],
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
econtrol1 <- trainControl(## 5-fold CV, repeat 5 times
  method = "repeatedcv",
  number = 5,
  ## repeated ten times
  repeats = 5,search = "random")
model_Keras_mlp<- caret::train(x=data1[,-1],y=(data1[,1]),
                                method=modelmlpkerasdropout,
                                metric = "MAE",## "Accuracy", "RMSE","Rsquared"
                                tuneLength = 11^3, ### 11 tunable parameters 11^2
                                # tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                                verbose=0,# verbose=1 is reporting the progress,o is sclience
                                trControl = econtrol1,importance = FALSE)





modelRSNNSmlpdecay <- list(label = "Multi-Layer Perceptron, multiple layers",
                                   library = "RSNNS",
                                   loop = NULL,
                                   type = c('Regression', 'Classification'),
                                   parameters = data.frame(parameter = c('layer1','layer2','layer3', 'decay'),
                                                           class = c('numeric','numeric','numeric', 'numeric'),
                                                           label = c('#Hidden Units layer1','#Hidden Units layer2','#Hidden Units layer3', 'Weight Decay')),
                                   grid = function(x, y, len = NULL, search = "grid"){
                                     if(search == "grid") {
                                       out <- expand.grid(layer1 = ((1:len) * 2) - 1, layer2 = 0, layer3 = 0,
                                                          decay = c(0, 10 ^ seq(-1, -4, length = len - 1)))
                                     } else {
                                       out <- data.frame(layer1 = sample(2:20, replace = TRUE, size = len),
                                                         layer2 = sample(c(0, 2:20), replace = TRUE, size = len),
                                                         layer3 = sample(c(0, 2:20), replace = TRUE, size = len),
                                                         decay = 10^runif(len, min = -5, max = 1))
                                     }
                                     out
                                   },
                                   fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                                     theDots <- list(...)
                                     theDots <- theDots[!(names(theDots) %in% c("size", "linOut"))]
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
                                     
                                     args <- list(x = x,
                                                  y = y,
                                                  learnFunc = "BackpropWeightDecay",
                                                  learnFuncParams = prms,
                                                  size = nodes,
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
                                   prob = function(modelFit, newdata, submodels = NULL) {
                                     out <- predict(modelFit, newdata)
                                     colnames(out) <- modelFit$obsLevels
                                     out
                                   },
                                   levels = function(x) x$obsLevels,
                                   tags = c("Neural Network","L2 Regularization"),
                                   sort = function(x) x[order(x$layer1, x$layer2, x$layer3, -x$decay),])


econtrol1 <- trainControl(## 5-fold CV, repeat 5 times
  method = "repeatedcv",
  number = 5,
  ## repeated ten times
  repeats = 5,search = "random")
model_RSNNS_mlp<- caret::train(x=data1[,-1],y=(data1[,1]),
                               method=modelRSNNSmlpdecay,
                               metric = "MAE",## "Accuracy", "RMSE","Rsquared"
                               tuneLength = 11^3, ### 11 tunable parameters 11^2
                               # tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                               verbose=0,# verbose=1 is reporting the progress,o is sclience
                               trControl = econtrol1,importance = FALSE)
