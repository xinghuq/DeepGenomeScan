\name{DLqvalues}
\alias{DLqvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DLqvalues(DL_data, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DL_data}{
%%     ~~Describe \code{DL_data} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
