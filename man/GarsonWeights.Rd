\name{GarsonWeights}
\alias{GarsonWeights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimating Garson weights from neural networks
%%  ~~function to do ... ~~
}
\description{Based on Gevrey et al. (2003), this function estimates the weights of the nodes/loci in the neural networks. For each hidden neuron i, multiply the absolute value of the hidden-output layer connection weight by the absolute value of the hidden-input layer connection weight. Do this for each input variable j. The following products Pij are obtained.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GarsonWeights(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The neural netwrok model from one of the h2o, RFCNN4, RSNNS models.
%%     ~~Describe \code{object} here~~
}
}
\details{Based on Gevrey et al. (2003), this function estimates the weights of the nodes in the neural networks. For each hidden neuron i, multiply the absolute value of the hidden-output layer connection weight by the absolute value of the hidden-input layer connection weight. Do this for each input variable j. The following products Pij are obtained.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ The nueral weights.
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
\keyword{Neural networks}% use one of  RShowDoc("KEYWORDS")

