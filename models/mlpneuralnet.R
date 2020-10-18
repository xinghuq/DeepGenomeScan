### now we use self-defined smooth RELU function
softplus <- function(x) log(1 + exp(x)) ### indeed this function learns well with old version of neuralnet package. return to the older version if necessary which will have less hyperparameter (difficult to overfit with older version)

### mlpneuralnet

## there is no output function option for the default version, but linear out choice determines the output function 
modelInfo<- list(label = "Neural Network",
                      library = "neuralnet",
                      loop = NULL,
                      type = c('Regression'),
                      parameters = data.frame(parameter = c('layer1', 'layer2', 'layer3',"activation1","linear.output"),
                                              class = c('numeric', 'numeric', 'numeric',"character","character"),
                                              label = c('Number of hidden Units in Layer 1', 'number of hidden Units in Layer 2', 'number of hidden Units in Layer 3',"Activation function in hidden layer","Activation function linear out choice")),
                      grid = function(x, y, len = NULL, search = "grid") {
                        afuncs=c("logistic", "tanh","softplus")
                        outputf=c("TRUE","FALSE")
                        if(search == "grid") {
                          out <- expand.grid(layer1 = ((1:len) * 2) - 1, layer2 = 0, layer3 = 0, activation1=c("logistic", "tanh","softplus"),linear.output=c("TRUE","FALSE"))
                        } else {
                          out <- data.frame(layer1 = sample(2:20, replace = TRUE, size = len),
                                            layer2 = sample(c(0, 2:20), replace = TRUE, size = len),
                                            layer3 = sample(c(0, 2:20), replace = TRUE, size = len),
                                            activation1=sample(afuncs, size = len, replace = TRUE),
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
                        
                        linear.output=c(param$linear.output)
                        neuralnet::neuralnet(form, algorithm="rprop+",data = dat,rep=1, hidden = nodes, stepmax = 1e8, learningrate.factor = list(minus = 0.5,plus = 1.2),act.fct=actf,linear.output=linear.output,...)
                      },
                      predict = function(modelFit, newdata, submodels = NULL) {
                        newdata <- newdata[, modelFit$model.list$variables, drop = FALSE]
                        neuralnet::compute(modelFit, covariate = newdata)$net.result[,1]
                      },
                 varImp = function(object, ...){
                   imps <- NeuralNetTools::olden(object,bar_plot =FALSE)
                   out <- data.frame(Overall = as.vector(imps))
                   rownames(out) <- names(imps)
                   out
                 },
                      prob = NULL,
                      tags = c("Neural Network"),
                      sort = function(x) x[order(x$layer1, x$layer2, x$layer3,x$activation1,x$linear.output),])
