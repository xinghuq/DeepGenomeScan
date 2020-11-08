\name{model_info}
\alias{model_info}
\docType{data}
\title{ Available models in DeepGenomeScan
%%   ~~ data name/kind ... ~~
}
\description{ The model list that shows the available model in DeepGenomeScan
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("model_info")}
\format{
  A data frame with 35 observations on the following 6 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{name}}{a factor with levels \code{Convolutional Neural Network} \code{Deep learning (Multi-Layer Perceptron, MLP)} \code{glmnet (Lasso and Elastic-Net Regularized Generalized Models)} \code{Linear Regression} \code{Logic Regression} \code{Model Averaged Neural Network} \code{Monotone Multi-Layer Perceptron Neural Network} \code{Multi-Layer Perceptron (MLP)} \code{Multi-Layer Perceptron, MLP, multiple layers} \code{Multilayer Perceptron Network by Stochastic Gradient Descent} \code{Multilayer Perceptron Network with Dropout} \code{Multilayer Perceptron Network with Weight Decay} \code{Multi-Layer Perceptron Neural Network} \code{Multi-Layer Perceptron, with multiple layers} \code{Multi-Step Adaptive MCP-Net} \code{Neural Network} \code{Neural Networks with Feature Extraction} \code{Principal Component Analysis} \code{Quantile Regression Neural Network} \code{Ridge Regression} \code{Ridge Regression with Variable Selection} \code{Stacked AutoEncoder Deep Neural Network} \code{Stochastic Gradient Boosting} \code{The Bayesian lasso} \code{The lasso} \code{Tree Models from Genetic Algorithms}}
    \item{\code{model}}{a factor with levels \code{avNNet} \code{blasso} \code{CNNsgd} \code{dnn} \code{evtree} \code{foba} \code{gbm} \code{glmnet} \code{glmnet_h2o} \code{lasso} \code{lm} \code{logreg} \code{mlp} \code{mlpFCNN4Rsgd} \code{mlph2o} \code{mlpKerasDecay} \code{mlpKerasDropout} \code{mlpML} \code{mlpneuralnet} \code{mlpneuralnet1} \code{mlpSGD} \code{mlpWeightDecay} \code{mlpWeightDecayML} \code{modelmlpkerasdropout} \code{modelRSNNSmlpdecay} \code{monmlp} \code{msaenet} \code{mxnet} \code{mxnetAdam} \code{neuralnet} \code{nnet} \code{pcaNNet} \code{pcr} \code{qrnn} \code{ridge}}
    \item{\code{type}}{a factor with levels \code{Classification, Regression} \code{Regression}}
    \item{\code{libraries}}{a factor with levels \code{} \code{deepnet} \code{elasticnet} \code{evtree} \code{FCNN4R} \code{FCNN4R, plyr} \code{foba} \code{gbm, plyr} \code{glmnet, Matrix} \code{h2o} \code{keras} \code{LogicReg} \code{monmlp} \code{monomvn} \code{msaenet} \code{mxnet} \code{neuralnet} \code{nnet} \code{pls} \code{qrnn} \code{RSNNS}}
    \item{\code{num_param}}{a factor with levels \code{alpha} \code{alpha, lambda} \code{alphas, nsteps, scale} \code{fraction} \code{hidden1, n.ensemble} \code{intercept} \code{k, lambda} \code{lambda} \code{layer1, layer2, layer3} \code{layer1, layer2, layer3,activation1,activation2,linear.output} \code{layer1, layer2, layer3,activation1,linear.output} \code{layer1, layer2, layer3, decay} \code{layer1,layer2,layer3, decay,activation1,activation2} \code{layer1, layer2, layer3, dropout, beta1, beta2, learningrate, activation} \code{layer1, layer2, layer3, hidden_dropout, visible_dropout} \code{layer1, layer2, layer3, learning.rate, momentum, dropout, activation} \code{ncomp} \code{nFilter,nStride,lambda, units1,units2,dropout, activation1,activation2,activation3} \code{n.hidden, penalty, bag} \code{n.trees, interaction.depth, shrinkage, n.minobsinnode} \code{size} \code{size, decay} \code{size, decay, bag} \code{size, dropout, batch_size, lr, rho, decay, activation} \code{size, l2reg, lambda, learn_rate, momentum, gamma, minibatchsz, repeats} \code{size, lambda, batch_size, lr, rho, decay, activation} \code{sparsity} \code{treesize, ntrees} \code{units1,units2, dropout1,dropout2,batch_size, lr, rho, decay, activation1,activation2,activation3} \code{units1,units2, l2reg, lambda, learn_rate,  momentum, gamma, minibatchsz, repeats,activation1,activation2} \code{units1,units2, l2reg, rho, activation}}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(model_info)
## maybe str(model_info) ; plot(model_info) ...
}
\keyword{datasets}
