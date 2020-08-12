##### CNN implementation using different CV methods and hyperparameters tuning. 
#####Calculation of variable importance for regression, all the pipeline are defined 

### define a cnn model used in caret to be able to implement the resampling process so that can estimate the performance of cnn, further extracts the feature importance
library(caret)### for ML calling functions and performance estimation
library(keras) ### for DL
library("tensorflow")
library("caretEnsemble")
### This will use python tensoflow, keras, ensure the tensorflow (2.0.0b1) is installed and is compatible with R version
###make sure the below work so that the configures are fine
backend(convert = TRUE)
keras::backend() ## or use tensorflow::backend()

## else you may need to clink python again to use keras or tensorflow as the backend
#library(reticulate)
##use_python("/usr/bin/python3")

### example of using one environmental factor to test CNN
data0=(cbind(simdata[[1]][,para[1]],simdata[[1]][,-(1:14)]))
# colnames(data0)[1]=para[[i]]
colnames(data0)[1]="enviri" ### This is more efficient than calling the name of the env factor when traing the model
data1=as.data.frame(apply(data0,2,normalize))

save.image(file = "CNN_test_data.RData")


nSNPs=ncol(data1[,-1]) 
nStride=3  # stride between convolutions
nFilter=32 # no. of convolutions

CNN_sgd <- list(label = "Convolutional Neural Network",
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
                varImp = NULL,
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


set.seed(123)
econtrol1 <- trainControl(## 5-fold CV, repeat 5 times
  method = "repeatedcv",
  number = 5,
  ## repeated ten times
  repeats = 5)
econtrol1 <- trainControl(## 5-fold CV, repeat 5 times
  method = "repeatedcv",
  number = 5,
  ## repeated ten times
  repeats = 5,search = "random")

CNNGrid <- expand.grid(nFilter=32,nStride=3,
                       lambda = c(0.001, 0.002, 0.003),
                       units1= c(2,4, 10),
                       units2= c(2,4, 10),
                       dropout = c(0.02,0.03),
                       activation1="relu",
                       activation2="relu",
                       activation3="linear")

### one value test
CNNGrid <- expand.grid(nFilter=32,nStride=3,
                       lambda = c(0.003),
                       units1= c(10),
                       units2= c(4),
                       dropout = c(0.02),
                       activation1="relu",
                       activation2="relu",
                       activation3="linear")

#### npote this is not a tenfor feeding to the CNN, this defined a keras object, though having the same dim with kerasR:: expend-dims (x1=keras::k_expand_dims(data1[,-1],axis=3))




################training CNN####################
#### result gives a model performance and the selected optimnal model#### with resamples stored 
model_Keras_CNN<- caret::train(x=data1[,-1],y=(data1[,1]),
                               method=CNN_sgd,
                               metric = "MAE",## "Accuracy", "RMSE","Rsquared"
                               #tuneLength = 91, ### 9 tunable parameters 9^2
                               tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                               verbose=0,# verbose=1 is reporting the progress,o is sclience
                               trControl = econtrol1,importance = FALSE)
model_Keras_CNN1<- caret::train(x=data1[,-1],y=(data1[,1]),
                                method=CNN_sgd,
                                metric = "MAE",## "Accuracy", "RMSE","Rsquared"
                                tuneLength = 91, ### 9 tunable parameters 9^2
                                # tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                                verbose=0,# verbose=1 is reporting the progress,o is sclience
                                trControl = econtrol1,importance = FALSE)

save(model_Keras_CNN1,file = "Test_env1_CNN_random_91.RData")
cnnImp <- varImp(model_Keras_CNN, scale = FALSE)##### estimating the variable importance using the resamples from the above training

write.csv(cnnImp$importance,file = paste0("sim_CNNkeras_env.csv"))

### the ablve is an example of CNNGrid search, this works
Sys.time()

#Gevrey, M., Dimopoulos, I., & Lek, S. (2003). Review and comparison of methods to study the contribution of variables in artificial neural network models. Ecological Modelling, 160(3), 249-264.

#Quinlan, J. (1992). Learning with continuous classes. Proceedings of the 5th Australian Joint Conference On Artificial Intelligence, 343-348.

###########################one should follow the below examnple to compile a working model for specific data, then replace the model in above pipleline###
################################# a typical example of CNN implementation ##################################
##########################################################working##############################################
#### define the CNN model
nSNP=ncol(simdata[[1]][,-(1:14)]) 
nStride=3  # stride between convolutions
nFilter=32 # no. of convolutions
x0=kerasR::expand_dims(data1[,-1],axis = 2) #[1]  640 1000    1
### construct CNN model, note the input shap should adopt the input tensor 
model_cnn_1 <- keras_model_sequential() %>%
  layer_conv_1d(filters = 32, kernel_size = c(3), strides=3,activation = "relu",
                input_shape = c(nSNP,1),kernel_regularizer = keras::regularizer_l2(0.03)) %>%
  layer_max_pooling_1d(pool_size = c( 2)) %>%
  layer_flatten() %>%
  layer_dropout(rate=0.5) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 10, activation = "linear")%>%
  layer_dense(units = 1)
#Compile the model
model_cnn_1 %>% compile(
  optimizer='sgd',
  loss='mean_squared_error',
  metrics='mean_squared_error')
summary(model_cnn_1)

model_cnn_1 %>% fit(
  x0,as.matrix(data1[,1]),
  epochs=10,
  batch_size=20)

#### predict and evaluation
pred<-keras:: predict_on_batch(model_cnn_1,x0)
performance_metrics=caret::postResample(pred, data1[,1])
model_in_R=keras::serialize_model(model_cnn_1) ### save in R

model_in_keras=keras::unserialize_model(model_in_R)## recover r into tensorflow model

model_cnn_1 %>% evaluate(x0, as.matrix(data1[,1]))

pred.value<-max.col(pred)-1

################################# a typical example of CNN implementation end##################################





#########copy backup


CNN_sgd <- list(label = "Convolutional Neural Network",
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
                                       lambda = c(0, 10 ^ seq(-1, -4, length = len - 1)),  units1 = ((1:len) * 2) - 1, units2 = ((1:len) * 1) - 1, dropout = seq(0, .7, length = len),
                                       activation1 = "relu",activation2="sigmoid",activation3="linear"
                    )
                  } else {
                    n <- nrow(x)
                    out <- data.frame(nFilter=sample(10:500, replace = TRUE, size = len),nStride=sample(2:20, replace = TRUE, size = len),
                                      lambda = 10^runif(len, min = -5, 1),units1 = sample(2:n*2, replace = TRUE, size = len),units2 = sample(1:as.integer(n/2)-1, replace = TRUE, size = len),
                                      dropout = runif(len, max = .7),
                                      activation = sample(
                                        afuncs, 
                                        size = len, 
                                        replace = TRUE
                                      )
                    )
                  }
                  out
                },
                
                fit = function(x, y, wts, last,lev,param, ...) {
                  require(dplyr)
                  K <- keras::backend()
                  K$clear_session()
                  # if(!is.matrix(x)) x <- as.matrix(x)
                  ### here should reshape the data
                  x=kerasR::expand_dims(x,axis=2)
                  
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
                varImp = NULL,
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
