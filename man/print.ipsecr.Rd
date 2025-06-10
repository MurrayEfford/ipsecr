\name{print.ipsecr}
\alias{print.ipsecr}
\alias{summary.ipsecr}
\alias{trim.ipsecr}

\title{ Print, Trim or Summarise ipsecr Object }

\description{

  Print results from fitting a spatially explicit capture--recapture model or
  generate a list of summary values.

}

\usage{
\method{print}{ipsecr} (x, newdata = NULL, alpha = 0.05, call = TRUE, ...)
\method{summary}{ipsecr} (object, newdata = NULL, alpha = 0.05, ...)
\method{trim}{ipsecr} (object, drop = c('call', 'proxyfn', 'mask', 'sim.lm'), keep = NULL)
}

\arguments{
  \item{x}{ \code{ipsecr} object output from \code{ipsecr.fit} }
  \item{object}{ \code{ipsecr} object output from \code{ipsecr.fit} }
  \item{newdata}{ optional dataframe of values at which to evaluate model }
  \item{alpha}{ alpha level }
  \item{call}{ logical; if TRUE the call is printed }
  \item{...}{ other arguments (not used) }
  \item{drop}{ character vector identifying components to be dropped } 
  \item{keep}{ character vector identifying components to be kept } 
}

\details{

Results from \code{print.ipsecr} are potentially complex and depend upon the analysis (see
below). Optional \code{newdata} should be a dataframe with a column for
each of the variables in the model. If \code{newdata} is missing then a
dataframe is constructed automatically.  Default \code{newdata} are for
a naive animal on the first occasion; numeric covariates are set to zero
and factor covariates to their base (first) level. Confidence intervals
are 100 (1 -- alpha) \% intervals.

\tabular{ll}{
call            \tab the function call (optional) \cr
version,time    \tab ipsecr version, date and time fitting started, and elapsed time \cr

Detector type   \tab `single', `multi', `proximity' etc. \cr
Detector number \tab number of detectors \cr
Average spacing \tab \cr
x-range         \tab \cr
y-range         \tab \cr
New detector type \tab as fitted when details$newdetector specified \cr
N animals       \tab number of distinct animals detected \cr
N detections    \tab number of detections \cr
N occasions     \tab number of sampling occasions \cr
Mask area       \tab \cr

Model          \tab model formula for each `real' parameter \cr
Fixed (real)   \tab fixed real parameters \cr
Detection fn   \tab detection function type (halfnormal or hazard-rate) \cr
Distribution    \tab spatial model (details$distribution) \cr
N parameters   \tab number of parameters estimated  \cr

Design points   \tab number of vertices and centre points \cr
Simulations per box \tab total number \cr
Beta parameters \tab coef of the fitted model, SE and confidence intervals \cr
vcov            \tab variance-covariance matrix of beta parameters  \cr
Real parameters \tab fitted (real) parameters evaluated at base levels of covariates \cr
}

}

\value{

The \code{summary} method constructs a list of outputs similar to those printed by the \code{print} method, but somewhat more concise and re-usable:

\tabular{ll}{
versiontime     \tab ipsecr version, and date and time fitting started \cr

traps           \tab detector summary \cr
capthist        \tab capthist summary \cr
mask            \tab mask summary \cr
modeldetails    \tab miscellaneous model characteristics \cr
coef            \tab table of fitted coefficients with CI \cr
predicted       \tab predicted values (`real' parameter estimates) \cr
}

}

\seealso{ 

\code{\link{ipsecr.fit}},
\code{\link[secr]{trim}}

}

\examples{

## load & print previously fitted null (constant parameter) model

print(ipsecrdemo)

summary(ipsecrdemo)

}

\keyword{ print }
