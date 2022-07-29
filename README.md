[![R build status](https://https://github.com/xinghuq/DeepGenomeScan//workflows/R-CMD-check/badge.svg)](https://github.com/xinghuq/DeepGenomeScan/)
[![Build Status](https://travis-ci.com/xinghuq/DeepGenomeScan.svg?branch=master)](https://travis-ci.com/xinghuq/DeepGenomeScan/)
[![Build status](https://ci.appveyor.com/api/projects/status/hqpwmgfeuerel48l?svg=true)](https://ci.appveyor.com/xinghuq/DeepGenomeScan/)
[![License:C](https://github.com/xinghuq/DeepGenomeScan/blob/master/License)](https://github.com/xinghuq/DeepGenomeScan/blob/master/License)


## DeepGenomeScan : A Deep Learning Approach for Whole Genome Scan (WGS) and Genome-wide Association Studies (GWAS)

  This package implements the genome scan and genome-wide association studies using deep neural networks (i.e, Multi-Layer Perceptron (MLP), Convolutional Neural Network (CNN)). DeepGenomeScan offers heuristic computational framework integrating different neural network architectures (i.e.,Multi-Layer Perceptron (MLP), convolutional neural network(CNN)) and robust resampling  methods, as well as the Model-Agnostic interpretation of feature importance for convolutional neural networks. DeepGenomeScan, in other words, deep learning for genome-wide scanning, is a deep learning approach for detecting signatures of natural selection and for performing various omics-based genome-wide association studies, such as GWAS, PWAS, TWAS, MWAS. The design makes the implemention user-friendly. It is compatible with most self-defined machine learning models (the self-defined models shuold be complete, including tunable parameters, fitted model, predicted model, examples can be found in our tutorial). Users can adopt the package's framework to study various ecological and evolutionary questions.

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

 install.packages("keras")
  
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
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
genotype_norm=as.data.frame(apply(genotype,2,normalize))
``````


### Setting the resampling method
``````{r}
econtrol1 <- trainControl(## 5-fold CV, repeat 5 times
  method = "adaptive_cv",
  number = 5,
  ## repeated ten times
  repeats = 5,search = "random")
set.seed(999)
options(warn=-1)
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
##### Calculating the p-values for SNPs and plot the SNP Manhattan plot

``````{r}

DLqvaluesarsine<-function(DL_data,K)
{
  loadings<-DL_data# [,1:as.numeric(K)]
  normdat <- apply(loadings, 2, normalize)
  asindat=apply(normdat,2, function(x) {asin(sqrt(x))})
  resmaha <- covRob(asindat, distance = TRUE, na.action= na.omit, estim="donostah")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_DL<-qval$qvalues
  padj <- p.adjust(reschi2test,method="bonferroni")
  return(data.frame(p.values=reschi2test, q.values=q.values_DL,padj=padj,mahaD=resmaha))
}



DLsim1=apply(out,2,normalize) #18

Simqvaluear=DLqvaluesarsine(DLsim1,1)
``````

## Manhattan plot

``````
ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), y=-log10(Simqvaluear[-which(Loci!="Neutral"),1])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), y=-log10(Simqvaluear[-which(Loci=="Neutral"),1]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("DNN -log10(p-value)") +ylim(c(0,100))+theme_bw()
plot(out$Overall,  ylab="SNP importance")

``````
### Package tutorial

More about how to construct your own model using DeepGenomeScan can be found in the [tutorial](https://xinghuq.github.io/DeepGenomeScan/index.html)

Welcome any [feedback](https://github.com/xinghuq/DeepGenomeScan/issues) and [pull request](https://github.com/xinghuq/DeepGenomeScan/pulls). 

## Version

The current version is 0.5.5 (Sep 2, 2020).

## Contact

qinxinghu@gmail.com

## Citation

Qin X, Chiang CWK, Gaggiotti OE. 2022. Deciphering signatures of natural selection via deep learning. Briefings in Bioinformatics. preprint (bioRxiv:2021.2005.2027.445973).
