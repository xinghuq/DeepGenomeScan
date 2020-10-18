\name{DeepGenomeScan}
\alias{DeepGenomeScan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Implementation of deep learning based genome scan
%%  ~~function to do ... ~~
}
\description{ Deep learning based genome scan with various ready-to-use models and self-defined models. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DeepGenomeScan(x, ...)
}
## Default S3 method:
DeepGenomeScan(
  genotype, 
  env,
  method = "modelRSNNSmlpdecay",
  preProcess = NULL,
  ...,
  weights = NULL,
  metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
  maximize = ifelse(metric %in% c("RMSE", "logLoss", "MAE"), FALSE, TRUE),
  trControl = trainControl(),
  tuneGrid = NULL,
  tuneLength = ifelse(trControl$method == "none", 1, 3),
seed = 123)

## S3 method for class 'formula'
DeepGenomeScan(form, data, ..., weights, subset, na.action = na.fail, contrasts = NULL,seed = 123)

## S3 method for class 'recipe'
DeepGenomeScan(
  genotype,
  data,
  method = "modelRSNNSmlpdecay",
  ...,
  metric = ifelse(is.factor(y_dat), "Accuracy", "RMSE"),
  maximize = ifelse(metric %in% c("RMSE", "logLoss", "MAE"), FALSE, TRUE),
  trControl = trainControl(),
  tuneGrid = NULL,
  tuneLength = ifelse(trControl$method == "none", 1, 3),
seed = 123)
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype}{The genotype matrix, or other omics-based variations matrix. This is a  matrix where samples are in rows and variations/features are in columns. This could be a simple matrix, data frame or other type (e.g. sparse matrix) but must have column names (see Details below). Preprocessing using the preProcess argument only supports matrices or data frames. When using the recipe method, _"genotype"_ should be an unprepared recipe object that describes the model terms (i.e. outcome, predictors, etc., e.g. the h2o model in our example) as well as any pre-processing that should be done to the data. This is an alternative approach to specifying the model. Note that, when using the recipe method, any arguments passed to preProcess will be ignored. See the links and example below for more details using recipes.
%%     ~~Describe \code{x} here~~
}
\item{env/traits}{The environmental variables/geographic coordinates, or the latent geographic genetic structure vairbales produced from KLFDAPC. For GWAS, this shoud be the traits or phenotypic variables. The data can be a numeric or factor vector corresponding to each sample.
%%     ~~Describe \code{x} here~~
}

 \item{\methods}{The read-to-use models or user defined models. A string specifying which model to use. A list of functions can also be passed for a custom model function, see our example tutorials.
%%     ~~Describe \code{\dots} here~~
}
  \item{\dots}{Arguments passed to the defined models (such as "mlpWeightDecayML"). Errors will occur if values for tuning parameters are passed here.
%%     ~~Describe \code{\dots} here~~
}



}
\details{
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
function (x, ...) 
{
    UseMethod("DeepGenomeScan")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
