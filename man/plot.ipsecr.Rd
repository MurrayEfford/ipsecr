\name{plot.ipsecr}
\alias{plot.ipsecr}

\title{ Plot Detection Functions }

\description{
  Plot detection functions using estimates of parameters in an ipsecr object.
}

\usage{

\method{plot}{ipsecr}(x, newdata = NULL, add = FALSE,
    sigmatick = FALSE, rgr = FALSE, limits = FALSE, alpha = 0.05,
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...)

}

\arguments{

  \item{x}{ an \code{ipsecr} object }
  \item{newdata}{ dataframe of data to form estimates }
  \item{add}{ logical to add curve(s) to an existing plot }
  \item{sigmatick}{ logical; if TRUE the scale parameter sigma is shown by a vertical line }
  \item{rgr}{ logical; if TRUE a scaled curve r.g(r) is plotted instead of g(r) }
  \item{limits}{ logical; if TRUE pointwise confidence limits are drawn }
  \item{alpha}{ alpha level for confidence intervals }
  \item{xval}{ vector of distances at for which detection to be plotted }
  \item{ylim}{ vector length 2 giving limits of y axis }
  \item{xlab}{ label for x axis }
  \item{ylab}{ label for y axis }
  \item{\dots}{ arguments to pass to \code{lines} }
  
}

\details{

\code{newdata} is usually NULL, in which case one curve is plotted for
each session and group. Otherwise, \code{predict.ipsecr} is used to form
estimates and plot a curve for each row in \code{newdata}.

If axis labels are not provided they default to `Distance (m)' and `Detection probability' or `Detection lambda'.

Approximate confidence limits for g(r) are calculated using a numerical
first-order delta-method approximation to the standard error at each
\code{xval}. The distribution of g(r) is assumed to be normal on the logit scale for non-hazard functions (detectfn 0:13). For hazard detection functions (detectfn 14:18) the hazard is assumed (from version 3.1.1) to be distributed normally on the log scale. Limits are back-transformed to the probability scale g(r).

}

\value{

\code{plot.ipsecr} invisibly returns a dataframe of the plotted values (or
a list of dataframes in the case that \code{newdata} has more than one
row).

}

\seealso{ \code{\link[secr]{Detection functions}}, \code{\link{plot}}, \code{\link{ipsecr}}, \code{\link[secr]{detectfnplot}} }

\examples{

plot (ipsecrdemo, xval = 0:100, ylim = c(0, 0.4))

}

\keyword{ hplot }

