\name{GarsonWeights_FCNN4R}
\alias{GarsonWeights_FCNN4R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimating Garson weights from FCNN4R models 
%%  ~~function to do ... ~~
}
\description{This function estimates the weights of loci based on Garson et al.(2013). 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GarsonWeights_FCNN4R(object, xnames = NULL, ynames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{FCNN4R model
%%     ~~Describe \code{object} here~~
}
  \item{xnames}{ optional, independent variable names
%%     ~~Describe \code{xnames} here~~
}
  \item{ynames}{optional, dependent variable names
%%     ~~Describe \code{ynames} here~~
}
}
\details{This function estimates the weights of neurons/loci based on Garson et al.(2013). Based on Gevrey et al. (2003), for each hidden neuron i, multiply the absolute value of the hidden-output layer connection weight by the absolute value of the hidden-input layer connection weight. Do this for each input variable j. The following products Pij are obtained.
%%  ~~ If necessary, more details than the description above ~~
}
\value{an array with the weights of each locus.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Gevrey, M., Dimopoulos, I., & Lek, S. (2003). Review and comparison of methods to study the contribution of variables in artificial neural network models. Ecological modelling, 160(3), 249-264.
%% ~put references to the literature/web site here ~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Neural networks }% use one of  RShowDoc("KEYWORDS")
