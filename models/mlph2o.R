### compile DL_h2o in DeepGenomeScan framework, for GUP version, please use H2O4GPU package
#### different libraries have different data format/processing 
##DL_h2o
modelInfo<- list(label = "DL_h2o",
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
                             out <- data.frame(units1 = sample(1:20, replace = TRUE, size = len),
                                               units2 = sample(1:20, replace = TRUE, size = len),
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
