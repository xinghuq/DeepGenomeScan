####modelRSNNSmlpdecay

modelInfo <- list(label = "Multi-Layer Perceptron, multiple layers",
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
                               out <- data.frame(layer1 = sample(2:20, replace = TRUE, size = len),
                                                 layer2 = sample(c(0, 2:20), replace = TRUE, size = len),
                                                 layer3 = sample(c(0, 2:20), replace = TRUE, size = len),
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
