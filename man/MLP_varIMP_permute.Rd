\name{MLP_varIMP_permute}
\alias{MLP_varIMP_permute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MLP variable importance using permutation 
%%  ~~function to do ... ~~
}
\description{MLP permutation variable importance. The permutation importance measured by the decrease in a model score (i.e., Mean Decrease Accuracy (MDA), Mean Decrease in RMSE) when a variable is randomly shuffled n times. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MLP_varIMP_permute(optmodel, feature_names = NULL, train_y = NULL, train_x = NULL, type = c("difference", "ratio"), nsim = 1, sample_size = NULL, sample_frac = NULL, verbose = FALSE, progress = "none", parallel = FALSE, paropts = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{optmodel}{ The optimal model used to estimate variable importance
%%     ~~Describe \code{optmodel} here~~
}
  \item{feature_names}{ The names of the variables
%%     ~~Describe \code{feature_names} here~~
}
  \item{train_y}{ The Y variable (dependent variable) used in regression
%%     ~~Describe \code{train_y} here~~
}
  \item{train_x}{The independent variable dataset
%%     ~~Describe \code{train_x} here~~
}
  \item{type}{Type of comparison "difference" or "ratio"
%%     ~~Describe \code{type} here~~
}
  \item{nsim}{ number of permutations
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
\details{In this implementation, the best model is determined and the orignal variable metrics are used as the baseline. Then the permutation variable performance metrics are tested using the best model as the training set. This procedure breaks the relationship between the variable and the target, thus the drop in the model score is indicative of how much the model depends on the variable.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Return a list of scores, including MLP model decrease in accuracy, the permutation metrics, and the baseline metrics.
%%  ~Describe the value returne
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
    baseline <- as.data.frame(t(caret::postResample(pred = keras::predict_on_batch(optmodel, 
        train_x), obs = train_y)))
    type <- match.arg(type)
    `\%compare\%` <- if (type == "difference") {
        `-`
    }
    else {
        `/`
    }
    permute_columns <- function(x, columns = NULL) {
        if (is.null(columns)) {
            stop("No columns specified for permutation.")
        }
        x[, columns] <- x[sample(nrow(x)), columns]
        x
    }
    sort_importance_scores <- function(x, decreasing) {
        x[order(x$Importance, decreasing = decreasing), ]
    }
    MLP_varIMP <- replicate(nsim, (plyr::llply(feature_names, 
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
            train_x_permuted <- permute_columns(train_x, columns = x)
            permuted <- as.data.frame(t(caret::postResample(pred = keras::predict_on_batch(optmodel, 
                train_x_permuted), obs = train_y)))
        })))
    varimp1 = apply(MLP_varIMP, MARGIN = 1, function(x) do.call(cbind, 
        x))
    varimp2 = do.call(rbind, varimp1)
    nm <- colnames(varimp2)
    RMSE = varimp2[, grepl("^RMSE", nm)]
    Rsquared = varimp2[, grepl("^Rsquared", nm)]
    MAE = varimp2[, grepl("^MAE", nm)]
    RMSE = as.data.frame(rowMeans(RMSE, na.rm = TRUE))
    Rsquared = as.data.frame(rowMeans(Rsquared, na.rm = TRUE))
    MAE = as.data.frame(rowMeans(MAE, na.rm = TRUE))
    rownames(RMSE) = feature_names
    rownames(Rsquared) = feature_names
    rownames(MAE) = feature_names
    MLP_SNPsIMP = cbind(RMSE, Rsquared, MAE)
    Decrease_acc = lapply(1:3, function(i) (baseline[, i] - MLP_SNPsIMP[, 
        i]))
    MLP_Decrease_acc = do.call(cbind, Decrease_acc)
    return(list(MLP_Decrease_acc = MLP_Decrease_acc, MLP_SNPsIMP = MLP_SNPsIMP, 
        baseline = baseline))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
