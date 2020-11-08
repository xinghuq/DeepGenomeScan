[![R build status](https://https://github.com/xinghuq/DeepGenomeScan//workflows/R-CMD-check/badge.svg)](https://github.com/xinghuq/DeepGenomeScan/)
[![Build Status](https://travis-ci.com/xinghuq/DeepGenomeScan.svg?branch=master)](https://travis-ci.com/xinghuq/DeepGenomeScan/)
[![Build status](https://ci.appveyor.com/api/projects/status/hqpwmgfeuerel48l?svg=true)](https://ci.appveyor.com/xinghuq/DeepGenomeScan/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## DeepGenomeScan : A Deep Learning Approach for Whole Genome Scan (WGS) and Genome-wide Association Studies (GWAS)

  This package implements the genome scan and genome-wide association studies using deep neural networks (i.e, Multi-Layer Perceptron (MLP), Convolutional Neural Network (CNN)). DeepGenomeScan offers heuristic computational framework integrating deep learning (i.e.,Multi-Layer Perceptron (MLP), convolutional neural network(CNN)), robust resampling and cross validations methods, as well as Model-Agnostic interpretation of feature importance for convolutional neural networks. DeepGenomeScan, in other words, deep learning for genome-wide scanning, is a deep learning approach for detecting variations under natural selection or omics-based association studies, such as GWAS, PWAS, TWAS, MWAS. The design makes the implemention user-friendly. It is compatible with most self-defined machine learning models (the self-defined models shuold be complete, including fit model, predictive model). Users can adopt the package's framework to study various ecological and evolutionary questions, not only constraining in biology field.  

## Install packages
`````{r}
library("devtools")

devtools::install_github("xinghuq/DeepGenomeScan")
devtools::install_github("xinghuq/CaretPlus/pkg/caret")

``````
## Dependencies and environment requirements

#### Note: Environment requirements: python should be installed and the python package of Keras and Tensorflow should also be installed and work properly with the system

## Checking the python environment 
``````{r}
library("rappdirs")
library("reticulate")
reticulate::use_python("/usr/bin/python3")
library(caret) ### for ML calling functions and performance estimation, users should use the modified version at xinghuq/CaretPlus/caret instead of the original version
library(keras)  
library("tensorflow")

checking if Tensorflow works properly
K0=keras::backend()

``````
``````{r}
requireNamespace("KLFDAPC")

 if (!requireNamespace("KLFDAPC", quietly=TRUE))

  devtools::install_github("xinghuq/KLFDAPC")
  
 if (!requireNamespace("DA", quietly=TRUE))
 
  devtools::install_github("xinghuq/DA")

if (!requireNamespace("keras", quietly=TRUE))

 Install.packages("keras")
  
 if (!requireNamespace("tensorflow", quietly=TRUE))
 
  install.packages("tensorflow")
 if (!requireNamespace("kerasR", quietly=TRUE))
 
  install.packages("kerasR")

``````

### library

```{r library,message = FALSE}
library(DeepGenomeScan)
library(caret)### for ML calling functions and performance estimation
library(keras) ### for DL
library("tensorflow")
library("caretEnsemble")
library(kerasR)
library("RSNNS")
library(NeuralNetTools)

```

### Example

### Preparing data
``````{r}
f <- system.file('extdata',package='DeepGenomeScan')
infile <- file.path(f, "sim1.csv")
sim_example=read.csv(infile)
genotype=sim_example[,-c(1:14)]
env=sim_example[,2:11]
str(sim_example)
``````


### Setting the resampling method
``````{r}
econtrol1 <- trainControl(## 5-fold CV, repeat 5 times
  method = "adaptive_cv",
  number = 5,
  ## repeated ten times
  repeats = 5,search = "random")
set.seed(999)
options(warn=-1
``````
### DeepGenomeScan with "mlph2o" model
```{r}

h2o_mlp<- DeepGenomeScan(as.matrix(genotype_norm),env$envir1,
                                  method="mlph2o",
                                  metric = "RMSE",## "Accuracy", "RMSE","Rsquared","MAE"
                                  tuneLength = 10, ### 11 tunable parameters 11^2
                                  # tuneGrid=CNNGrid, ### or search 100 combinations of parameters using random tuneLength=100
                                  trControl = econtrol1)

#### There is a model specific varIMP for this model
varImp(h2o_mlp,scale = FALSE)

out <- as.data.frame(h2o::h2o.varimp(h2o_mlp$finalModel), stringsAsFactors = TRUE)
colnames(out)[colnames(out) == "scaled_importance"] <- "Overall"
rownames(out) <- out$names

``````
##### Plot the SNP importance scores
``````{r}
plot(out$Overall,  ylab="SNP importance")

``````
Welcome any [feedback](https://github.com/xinghuq/DeepGenomeScan/issues) and [pull request](https://github.com/xinghuq/DeepGenomeScan/pulls). 

## Version

The current version is 0.5.5 (Sep 2, 2020).

## Contact

qinxinghu@gmail.com

## Citation

Qin, X. 2020. DeepGenomeScan. v0.5.5.

Qin, X., Gaggiotti, EO., Chiang, C. 2020. Detecting natural selection via deep learning. In submission.
