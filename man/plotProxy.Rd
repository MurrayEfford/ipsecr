\name{plotProxy}
\alias{plotProxy}

\title{

Simulate and plot the relationship between a parameter and its designated proxy

}

\description{

As described in the vignette, each parameter is matched to a proxy value computed 
cheaply from the rawdata by the proxy function. This function provides a visual 
check on that relationship.

}

\usage{

plotProxy(parameter = "sigma", proxyfn = proxy.ms, traps, mask, detectfn = "HHN", 
    noccasions = 5, basepar = list(), xvals = NULL, nrepl = 100, add = FALSE, 
    trend = TRUE, points = FALSE, boxplot = TRUE, boxplotargs = list(), 
    link = "log", details = NULL, ...)
    
}

\arguments{
  \item{parameter}{character parameter of interest}
  \item{proxyfn}{function to compute vector of proxy values from a capthist object}
  \item{traps}{traps object}
  \item{mask}{habitat mask object}
  \item{detectfn}{numeric or character code for detection function (see \link{detectfn})}
  \item{noccasions}{integer number of sampling occasions}
  \item{basepar}{named list with central values of parameters}
  \item{xvals}{specified values of focal paramater to simulate (optional)}
  \item{nrepl}{integer number of simulations}
  \item{add}{logical; if TRUE any plot is added to an existing plot}
  \item{trend}{logical; if TRUE a least-squares trend line is plotted}
  \item{points}{logical; if TRUE each simulated value is plotted}
  \item{boxplot}{logical; if TRUE a boxplots is plotted for each level of the focal parameter}
  \item{boxplotargs}{list of arguments for \code{\link{boxplot}} (optional)}
  \item{link}{character link function for transforming parameter x-axis}
  \item{details}{not used}
  \item{\dots}{other arguments passed to plot()}
}

\details{

Simulation uses the internal functions \code{\link{simpop}} and \code{\link{simCH}}.

\code{boxplotargs} may be used to override or add to the arguments of \code{\link{boxplot}}.

This version of \code{plotProxy()} does not allow for interference (NT) and assumes a simple SECR model with only 3 or 4 coefficients corresponding to density D and the parameters of the detection model (lambda0 or g0, sigma and possibly z). 

Matching of proxies at the level of `beta' parameters may be enabled in a future version.

}

\value{

The simulated proxy values are returned invisibly as a matrix (nrepl x nlevels).

}

\seealso{

\code{\link{proxy.ms}}

}

\examples{

trps <- traps(captdata)
msk <- make.mask(trps, buffer = 100)
base <- list(D = 5, lambda0 = 0.2, sigma = 25)
out <- plotProxy (parameter = 'D', traps = trps, mask = msk,
     basepar = base, boxplotargs = list(col='orange'), nrepl = 50)

}
