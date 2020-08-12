\name{CNN_sgd}
\alias{CNN_sgd}
\docType{data}
\title{ Convolutional Neural Network (CNN) compiling object
%%   ~~ data name/kind ... ~~
}
\description{ This is an model compiling object used as an example for constructing CNN with CaretPlus in DeepGenomeScan framework
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("CNN_sgd")}
\format{
  The format is:
List of 14
 $ label     : chr "Convolutional Neural Network"
 $ library   : chr "keras"
 $ loop      : NULL
 $ type      : chr [1:2] "Regression" "Classification"
 $ parameters:'data.frame':	9 obs. of  3 variables:
  ..$ parameter: Factor w/ 9 levels "activation1",..: 6 7 5 8 9 4 1 2 3
  ..$ class    : Factor w/ 2 levels "character","numeric": 2 2 2 2 2 2 1 1 1
  ..$ label    : Factor w/ 9 levels "Activation Function convol layer",..: 8 9 7 5 6 4 1 2 3
 $ grid      :function (x, y, len = NULL, search = "grid")  
  ..- attr(*, "srcref")= 'srcref' int [1:8] 13 24 27 17 24 17 13 27
  .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x50b9ce8> 
 $ fit       :function (x, y, wts, last, lev, param, ...)  
  ..- attr(*, "srcref")= 'srcref' int [1:8] 29 23 101 17 23 17 29 101
  .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x50b9ce8> 
 $ predict   :function (modelFit, newdata, submodels = NULL)  
  ..- attr(*, "srcref")= 'srcref' int [1:8] 102 27 116 17 27 17 102 116
  .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x50b9ce8> 
 $ prob      :function (modelFit, newdata, submodels = NULL)  
  ..- attr(*, "srcref")= 'srcref' int [1:8] 117 25 126 17 25 17 117 126
  .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x50b9ce8> 
 $ varImp    : NULL
 $ tags      : chr [1:2] " CNN" "L2 Regularization"
 $ sort      :function (x)  
  ..- attr(*, "srcref")= 'srcref' int [1:8] 129 24 129 66 24 66 129 129
  .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x50b9ce8> 
 $ notes     : chr "After `train` completes, the keras model object is serialized so that it can be used between R session. When pr"| __truncated__
 $ check     :function (pkg)  
  ..- attr(*, "srcref")= 'srcref' int [1:8] 141 25 150 17 25 17 141 150
  .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x50b9ce8> 
}
\details{ CNN model compiling object
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(CNN_sgd)
## maybe str(CNN_sgd) ; plot(CNN_sgd) ...
}
\keyword{datasets}
