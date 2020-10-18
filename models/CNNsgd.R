
###CNNsgd

modelInfo <- list(label = "Convolutional Neural Network",
                library = "keras",
                loop = NULL,
                type = c('Regression', "Classification"),
                parameters = data.frame(
                  parameter = c("nFilter","nStride",'lambda', "units1","units2","dropout",
                                "activation1","activation2","activation3"),
                  class = c(rep('numeric', 6), "character","character","character"),
                  label = c("no. of convolutions","stride between convolutions",'L2 Regularization', "hidden_neurons1","hidden_neurons2",
                            "drop_out_rates","Activation Function convol layer","Activation function dense hidden layer","Activation function layer to output")
                ),
                
                grid = function(x, y, len = NULL, search = "grid") {
                  afuncs <- c("sigmoid", "relu", "tanh","linear")
                  if(search == "grid") {
                    out <- expand.grid(nFilter=as.integer(nrow(x)/10),nStride=3,
                                       lambda = c(0, 10 ^ seq(-1, -4, length = len - 1)),  units1 = ((1:len) * 2) - 1, units2 = ((1:len) * 1) - 1, dropout = seq(0, .5, length = len),
                                       activation1 = "relu",activation2="sigmoid",activation3="linear"
                    )
                  } else {
                    n <- nrow(x)
                    out <- data.frame(nFilter=sample(10:500, replace = TRUE, size = len),nStride=sample(2:20, replace = TRUE, size = len),
                                      lambda = 10^runif(len, min = -5, 1),units1 = sample(2:n*2, replace = TRUE, size = len),units2 = sample(1:as.integer(n/2)-1, replace = TRUE, size = len),
                                      dropout = runif(len, max = .5),activation1 = sample(afuncs,  size = len, replace = TRUE),activation2 = sample(afuncs,  size = len, replace = TRUE), activation3 = sample(afuncs,  size = len, replace = TRUE))
                  }
                  out
                },
                
                fit = function(x, y, wts, last,lev,param, ...) {
                  require(dplyr)
                  K <- keras::backend()
                  K$clear_session()
                  # if(!is.matrix(x)) x <- as.matrix(x)
                  ########## ### here should reshape the data################################# Note: the key step feed to the cnn
                  x=kerasR::expand_dims(x,axis=2)
                  ################  ###### also this define the shape of the tensor######  argument: input_shape    input_shape = c(nSNP,1)          
                  nSNP=dim(x)[2]
                  model<-keras::keras_model_sequential() %>%
                    keras::layer_conv_1d(filters = param$nFilter, kernel_size = c(3), strides=param$nStride,activation =as.character(param$activation1),
                                         input_shape = c(nSNP,1),kernel_regularizer = keras::regularizer_l2(param$lambda)) %>%
                    layer_max_pooling_1d(pool_size = c( 2)) %>% # add pooling layer: takes maximum of two consecutive values
                    layer_flatten() %>%
                    layer_dropout(rate=param$dropout) %>%
                    layer_dense(units = param$units1, activation =as.character(param$activation2)) %>%
                    layer_dense(units = param$units2, activation =as.character(param$activation3))
                  
                  
                  #                    model<-keras::keras_model_sequential() %>%
                  #                      layer_conv_1d(filters = param$nFilter, kernel_size = c(3), strides=param$nStride,activation =as.character(param$activation1),
                  #                                    input_shape = c(1000,1),kernel_regularizer = keras::regularizer_l2(param$lambda)) %>%
                  #                      layer_max_pooling_1d(pool_size = c( 2)) %>% # add pooling layer: takes maximum of two consecutive values
                  #                      layer_flatten() %>%
                  #                      layer_dropout(rate=param$dropout) %>%
                  #                      layer_dense(units = param$units1, activation =as.character(param$activation2)) %>%
                  #                      layer_dense(units = param$units2, activation =as.character(param$activation3))%>% #### activation function of last layer
                  #                      layer_dense(units = 1) ### this should be the same as the input shape: input_shape = c(nSNP,1)
                  #more complex model
                  #  model<-keras::keras_model_sequential() %>%
                  #    keras::layer_conv_1d(filters=128, kernel_size = 3, strides=3,activation =as.character(param$activation1),
                  #                  input_shape = c(nSNPs, 1),kernel_regularizer = keras::regularizer_l2(param$lambda)) %>%
                  #    layer_max_pooling_1d(pool_size = 2) %>%
                  #    layer_conv_1d(filters = 64, kernel_size =3, activation = as.character(param$activation1)) %>%
                  #    layer_flatten() %>%
                  #    layer_dropout(rate=param$dropout) %>%
                  #    layer_dense(units =param$units, activation = as.character(param$activation2))
                  #Compile the model. Note that this linked above %>% 
                  #   model%>% compile(
                  #     optimizer='sgd',
                  #     loss='mean_squared_error',
                  #     metrics='mean_squared_error')
                  
                  
                  if(is.factor(y)) {
                    y <- class2ind(y)
                    model %>% 
                      layer_dense(units = 1) %>%  #### activation function of last layer
                      keras::compile(
                        loss = "categorical_crossentropy",
                        optimizer = "sgd",
                        metrics = "accuracy"
                      )
                  } else {
                    model %>%
                      layer_dense(units = 1) %>% #### activation function of last layer
                      keras::compile(
                        loss = "mean_squared_error",
                        optimizer = "sgd",
                        metrics = "mean_squared_error"
                      )
                  }
                  
                  model %>% keras::fit(
                    x = x, 
                    y = y,
                    kernel_regularizer = keras::regularizer_l2(param$lambda),
                    ...
                  )
                  if(last)
                    model <- keras::serialize_model(model)
                  list(object = model)
                },
                predict = function(modelFit, newdata, submodels = NULL) {
                  if(inherits(modelFit, "raw"))
                    modelFit$object <- keras::unserialize_model(modelFit)
                  # if(!is.matrix(newdata)) 
                  # newdata <- as.matrix(newdata)
                  newdata=kerasR::expand_dims(newdata,axis=2) 
                  out <- keras::predict_on_batch(modelFit$object, newdata)
                  ## check for model type
                  if(ncol(out) == 1) {
                    out <- out[, 1]
                  } else {
                    out <- modelFit$obsLevels[apply(out, 1, which.max)]
                  }
                  out
                },
                prob =  function(modelFit, newdata, submodels = NULL) {
                  if(inherits(modelFit, "raw"))
                    modelFit$object <- keras::unserialize_model(modelFit)
                  # if(!is.matrix(newdata)) 
                  #newdata <- as.matrix(newdata)
                  newdata=kerasR::expand_dims(newdata,axis=2) 
                  out <- keras::predict_on_batch(modelFit$object, newdata)
                  colnames(out) <- modelFit$obsLevels
                  as.data.frame(out, stringsAsFactors = TRUE)
                },
                varImp = NULL,### CNN importance will estimate separately
                tags = c(" CNN", "L2 Regularization"),
                sort = function(x) x[order(x$nFilter, -x$lambda),],
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