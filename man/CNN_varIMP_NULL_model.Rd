\name{CNN_varIMP_NULL_model}
\alias{CNN_varIMP_NULL_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CNN variable importance using NULL model
%%  ~~function to do ... ~~
}
\description{CNN NULL variable importance. The NULL importance measured by the decrease in a model score (i.e., Mean Decrease Accuracy (MDA), Mean Decrease in RMSE) when a variable is set as NULL variable.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CNN_varIMP_NULL_model(optmodel, feature_names = NULL, train_y = NULL, train_x = NULL, smaller_is_better = NULL, type = c("difference", "ratio"), nsim = 1, sample_size = NULL, sample_frac = NULL, verbose = FALSE, progress = "none", parallel = FALSE, paropts = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{optmodel}{The optimal model
%%     ~~Describe \code{optmodel} here~~
}
  \item{feature_names}{The variable names, should be the same as the dim of X
%%     ~~Describe \code{feature_names} here~~
}
  \item{train_y}{ The Y variable (dependent variable) used in regression
%%     ~~Describe \code{train_y} here~~
}
  \item{train_x}{The independent variable dataset
%%     ~~Describe \code{train_x} here~~
}
  \item{smaller_is_better}{
%%     ~~Describe \code{smaller_is_better} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{nsim}{
%%     ~~Describe \code{nsim} here~~
}
  \item{sample_size}{
%%     ~~Describe \code{sample_size} here~~
}
  \item{sample_frac}{
%%     ~~Describe \code{sample_frac} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{progress}{
%%     ~~Describe \code{progress} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{paropts}{
%%     ~~Describe \code{paropts} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{This function estimates the similar model agnostic feature importance but rather using time consuming permutation importance, it uses the Null model importance. Likewise, the best model is determined and the orignal variable metrics are used as the baseline. Then the NULL variable performance metrics are tested using the best model as the training set. This procedure directly set the target variable as NULL variable (Zero variance), thus the drop in the model score is indicative of how much the model depends on the variable.

%%  ~~ If necessary, more details than the description above ~~
}
\value{Return a list of scores, including CNN model decrease in accuracy, the permutation metrics, and the baseline metrics.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
function (optmodel, feature_names = NULL, train_y = NULL, train_x = NULL, 
    smaller_is_better = NULL, type = c("difference", "ratio"), 
    nsim = 1, sample_size = NULL, sample_frac = NULL, verbose = FALSE, 
    progress = "none", parallel = FALSE, paropts = NULL, ...) 
{
    require(caret)
    x_cnn = kerasR::expand_dims(train_x, axis = 2)
    baseline <- as.data.frame(t(caret::postResample(pred = keras::predict_on_batch(optmodel, 
        x_cnn), obs = train_y)))
    type <- match.arg(type)
    `\%compare\%` <- if (type == "difference") {
        `-`
    }
    else {
        `/`
    }
    NULL_columns <- function(x, columns = NULL) {
        if (is.null(columns)) {
            stop("No columns specified for permutation.")
        }
        x[, columns] <- rep(1, nrow(x))
        x
    }
    sort_importance_scores <- function(x, decreasing) {
        x[order(x$Importance, decreasing = decreasing), ]
    }
    CNN_varIMP <- replicate(nsim, (plyr::llply(feature_names, 
        .progress = "none", .parallel = parallel, .paropts = paropts, 
        .fun = function(x) {
            if (verbose && !parallel) {
                message("Computing variable importance for ", 
                  x, "...")
            }
            if (!is.null(sample_size)) {
                ids <- sample(length(train_y), size = sample_size, 
                  replace = FALSE)
                train_x <- train_x[ids, ]
                train_y <- train_y[ids]
            }
            train_x_NULL <- NULL_columns(train_x, columns = x)
            x_tensor = kerasR::expand_dims(train_x_NULL, axis = 2)
            permuted <- as.data.frame(t(caret::postResample(pred = keras::predict_on_batch(optmodel, 
                x_tensor), obs = train_y)))
        })))
    CNN_SNPsIMP = do.call(rbind, CNN_varIMP)
    rownames(CNN_SNPsIMP) = feature_names
    decrease_acc = lapply(1:3, function(i) baseline[, i] - CNN_SNPsIMP[, 
        i])
    CNN_Decrease_acc = do.call(cbind, decrease_acc)
    return(list(CNN_Decrease_acc = CNN_Decrease_acc, CNN_SNPsIMP = CNN_SNPsIMP, 
        baseline = baseline))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
