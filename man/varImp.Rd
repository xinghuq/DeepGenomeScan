\name{varImp}
\alias{varImp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variable importance for DeepGenomeScan models
%%  ~~function to do ... ~~
}
\description{Variable importance for neural networks in DeepGenomeScan.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## S3 method for class 'neuralnet'
varImp(object, ...)
## S3 method for class 'h2o'
varImp(object, ...)
## S3 method for class 'kears'
varImp(object, ...)
## S3 method for class 'FCNN4R'
varImp(object, ...)
## S3 method for class 'RSNNS'
varImp(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The final fitted model or the model trained from DeepGenomeScan
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{parameters to pass to the specific varImp methods
%%     ~~Describe \code{\dots} here~~
}
}
\details{The varible importance estimated here based on various methods. In terms of MLP from "neuralnet", "RSNNS", we used Olden et al (2003), which uses cas the sum of the product of raw input-hidden, hidden-output connection weights, proposed by Olden et al. 2004.  In terms of the CNN and MLP importance extimated from Keras and Tensorflow libearies, we used themodel agnostic feature importance based on Ribeiro et al, 2018. Specifical CNN and MLP importance can see function  CNN_varIMP_NULL_model, MLP_varIMP_NULL_model.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Beck, M.W. 2018. NeuralNetTools: Visualization and Analysis Tools for Neural Networks. Journal of Statistical Software. 85(11):1-20.

Garson, G.D. 1991. Interpreting neural network connection weights. Artificial Intelligence Expert. 6(4):46-51.

Goh, A.T.C. 1995. Back-propagation neural networks for modeling complex systems. Artificial Intelligence in Engineering. 9(3):143-151.

Olden, J.D., Jackson, D.A. 2002. Illuminating the 'black-box': a randomization approach for understanding variable contributions in artificial neural networks. Ecological Modelling. 154:135-150.

Olden, J.D., Joy, M.K., Death, R.G. 2004. An accurate comparison of methods for quantifying variable importance in artificial neural networks using simulated data. Ecological Modelling. 178:389-397.

Ribeiro, Marco Tulio, Sameer Singh, and Carlos Guestrin. "Model-agnostic interpretability of machine learning." arXiv preprint arXiv:1606.05386 (2016).
Fisher, Aaron, Cynthia Rudin, and Francesca Dominici. “Model Class Reliance: Variable importance measures for any machine learning model class, from the ‘Rashomon’ perspective.” http://arxiv.org/abs/1801.01489 (2018).
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, ...) 
{
    UseMethod("varImp")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
