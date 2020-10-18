\name{DLqvalues}
\alias{DLqvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculating q-values from importance values
%%  ~~function to do ... ~~
}
\description{ The 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DLqvalues(DL_data, K, estimation="auto")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DL_data}{The variable importance. Row names are variable names, columns are the importance obtained from indenpendent factors
%%     ~~Describe \code{DL_data} here~~
}
  \item{K}{ The degrees of freedom (non-negative, but can be non-integer) used in estimating the q-values.
%%     ~~Describe \code{K} here~~
}
\item{estimation}{The estimation method. A character string specifying the robust estimator to be used. The choices are: "mcd" for the Fast MCD algorithm of Rousseeuw and Van Driessen, "weighted" for the Reweighted MCD, "donostah" for the Donoho-Stahel projection based estimator, "M" for the constrained M estimator provided by Rocke, "pairwiseQC" for the orthogonalized quadrant correlation pairwise estimator, and "pairwiseGK" for the Orthogonalized Gnanadesikan-Kettenring pairwise estimator. The default "auto" selects from "donostah", "mcd", and "pairwiseQC" with the goal of producing a good estimate in a reasonable amount of time.}
}
\details{The covRob function selects a robust covariance estimator that is likely to provide a good estimate in a reasonable amount of time. Presently this selection is based on the problem size. The Donoho-Stahel estimator is used if there are less than 1000 observations and less than 10 variables or less than 5000 observations and less than 5 variables. If there are less than 50000 observations and less than 20 variables then the MCD is used. For larger problems, the Orthogonalized Quadrant Correlation estimator is used.

The MCD and Reweighted-MCD estimates (estim = "mcd" and estim = "weighted" respectively) are computed using the covMcd function in the robustbase package. By default, covMcd returns the reweighted estimate; the actual MCD estimate is contained in the components of the output list prefixed with raw.

The M estimate (estim = "M") is computed using the covMest function in the rrcov package. For historical reasons the Robust Library uses the MCD to compute the initial estimate.

The Donoho-Stahel (estim = "donostah") estimator is computed using the CovSde function provided in the rrcov package.

The pairwise estimators (estim = "pairwisegk" and estim = "pairwiseqc") are computed using the CovOgk function in the rrcov package.

Version 0.3-8 of the Robust Library: all of the functions origianlly contributed by the S-Plus Robust Library have been replaced by dependencies on the robustbase and rrcov packages. Computed results may differ from earlier versions of the Robust Library. In particular, the MCD estimators are now adjusted by a small sample size correction factor. Additionally, a bug was fixed where the final MCD covariance estimate produced with estim = "mcd" was not rescaled for consistency.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Results return p-values and q-values for variables
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
R. A. Maronna and V. J. Yohai (1995) The Behavior of the Stahel-Donoho Robust Multivariate Estimator. Journal of the American Statistical Association 90 (429), 330–341.

P. J. Rousseeuw and K. van Driessen (1999) A fast algorithm for the minimum covariance determinant estimator. Technometrics 41, 212–223.

D. L. Woodruff and D. M. Rocke (1994) Computable robust estimation of multivariate location and shape on high dimension using compound estimators. Journal of the American Statistical Association, 89, 888–896.

R. A. Maronna and R. H. Zamar (2002) Robust estimates of location and dispersion of high-dimensional datasets. Technometrics 44 (4), 307–317.

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
function (DL_data, K) 
{
    loadings <- DL_data
    resscale <- apply(loadings, 2, scale)
    resmaha <- covRob(resscale, distance = TRUE, na.action = na.omit, 
        estim = "donostah")$dist
    lambda <- median(resmaha)/qchisq(0.5, df = K)
    reschi2test <- pchisq(resmaha/lambda, K, lower.tail = FALSE)
    qval <- qvalue(reschi2test)
    q.values_DL <- qval$qvalues
    padj <- p.adjust(reschi2test, method = "bonferroni")
    return(data.frame(p.values = reschi2test, q.values = q.values_DL, 
        padj = padj))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
