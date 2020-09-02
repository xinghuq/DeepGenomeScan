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

  devtools::Install_github("xinghuq/KLFDAPC")
  
 if (!requireNamespace("DA", quietly=TRUE))
 
  devtools::install_github("xinghuq/DA")

if (!requireNamespace("keras", quietly=TRUE))

 Install.packages("keras")
  
 if (!requireNamespace("tensorflow", quietly=TRUE))
 
  install.packages("tensorflow")
 if (!requireNamespace("kerasR", quietly=TRUE))
 
  install.packages("kerasR")

``````
### Example

``````{r}
f <- system.file('extdata',package='DeepGenomeScan')

infile <- file.path(f, "Test_env1_CNN_random_91.RData")
``````
##Note the labels below is random
``````{r}
y1=rep(1,times=1736)
y2=rep(2,times=2000)

y=rbind(as.matrix(y1),as.matrix(y2))

y=as.factor(y)

``````
``````{r}
# Using gaussan kernel
# This will take longer than PCA, denpending on the number of samples and n.pcs. We will not show the results here. Users can test on their own clusters

virus_klfdapc=KLFDAPC(infile,y,kernel=kernlab::rbfdot(sigma = 0.5),r=3,snp.id=NULL, maf=0.05, missing.rate=0.05,n.pc=10,tol=1e-30, num.thread=2,metric = "plain",prior = NULL)

showfile.gds(closeall=TRUE)
``````
##### Plot the reduced features
``````{r}
plot(virus_klfdapc$KLFDAPC$Z, virus_klfdapc$KLFDAPC$Z, col=as.integer(y), xlab="KLFDA 2", ylab="KLFDA 1")
legend("bottomright", legend=levels(y), pch="o", col=1:nlevels(y))

``````
Welcome any [feedback](https://github.com/xinghuq/DeepGenomeScan/issues) and [pull request](https://github.com/xinghuq/DeepGenomeScan/pulls). 

## Version

The current version is 0.5.5 (Sep 2, 2020).

## Citation

Qin. X. 2020. DeepGenomeScan. v0.5.5.
