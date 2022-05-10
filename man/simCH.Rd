\name{simCH}
\alias{simCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simCH(traps, popn, detectfn, detectpar, noccasions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traps}{
%%     ~~Describe \code{traps} here~~
}
  \item{popn}{
%%     ~~Describe \code{popn} here~~
}
  \item{detectfn}{
%%     ~~Describe \code{detectfn} here~~
}
  \item{detectpar}{
%%     ~~Describe \code{detectpar} here~~
}
  \item{noccasions}{
%%     ~~Describe \code{noccasions} here~~
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
function (traps, popn, detectfn, detectpar, noccasions) 
{
    K <- nrow(traps)
    usge <- usage(traps)
    if (is.null(usge)) {
        usge <- matrix(1, K, noccasions)
    }
    detectcode <- switch(detector(traps)[1], single = -1, multi = 0, 
        proximity = 1, count = 2, capped = 8, 9)
    if (detectcode == 9) 
        stop("unsupported detector type")
    lambdak <- detectpar[["lambdak"]][1]
    if (is.null(lambdak)) 
        lambdak <- -1
    temp <- CHcpp(as.matrix(popn), as.matrix(traps), as.matrix(usge), 
        as.integer(detectfn), as.integer(detectcode), unlist(detectpar[parnames(detectfn)]), 
        as.double(lambdak), 0, 0, 0)
    if (temp$resultcode != 0) {
        stop("simulated detection failed, code ", temp$resultcode)
    }
    npop <- nrow(popn)
    w <- array(temp$value, dim = c(noccasions, K, npop), dimnames = list(1:noccasions, 
        NULL, 1:npop))
    w <- aperm(w, c(3, 1, 2))
    if (modelnontarget) {
        nontarget <- temp$nontarget
    }
    else {
        nontarget <- NULL
    }
    w <- w[apply(w, 1, sum) > 0, , , drop = FALSE]
    class(w) <- "capthist"
    attr(w, "nontarget") <- nontarget
    traps(w) <- traps
    w
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
