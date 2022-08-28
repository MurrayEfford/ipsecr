\name{ipsecr-package}
\alias{ipsecr-package}
\alias{ipsecr}

\docType{package}

\title{Spatially Explicit Capture--Recapture by Inverse Prediction }

\description{

  Functions to estimate the density and size of a spatially
  distributed animal population sampled with an array of passive
  detectors, such as traps. \pkg{ipsecr} addresses `difficult' models that 
  strictly cannot be fitted by maximum likelihood in the larger parent package 
  \pkg{secr} (Efford 2022). The classic example concerns data from single-catch traps.

}

\details{

\tabular{ll}{
Package: \tab ipsecr\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2022-08-28\cr
License: \tab GNU General Public License Version 2 or later\cr
}

Spatially explicit capture--recapture is a set of methods for studying
marked animals distributed in space. Data comprise the locations of
detectors (described in an object of class `traps'), and the 
detection histories of individually marked
animals. Individual histories are stored in an object of class
`capthist' that includes the relevant `traps' object.

Models for population density (animals per hectare) and detection are
defined in \pkg{ipsecr} using symbolic formula notation (but only constant models are allowed for detection parameters in \pkg{ipsecr} 1.1) Density models
may include spatial trend. Habitat is distinguished from
nonhabitat with an object of class `mask'.

Models are fitted in \pkg{ipsecr} by simulation and inverse prediction (Efford 2004).
A model fitted with \code{\link{ipsecr.fit}} is an object
of class \code{ipsecr}. Generic methods (plot, print, summary, etc.) are
provided.

A link at the bottom of each help page takes you to the
help index. Several vignettes complement the help pages:

The analyses in \pkg{ipsecr} extend those available in the software
Density (see
\href{https://www.otago.ac.nz/density/}{www.otago.ac.nz/density/}
for the most recent version of Density).  Help is available on the
`DENSITY | secr' forum at
\href{http://www.phidot.org/forum/}{www.phidot.org} and
the Google group \href{https://groups.google.com/forum/#!forum/secrgroup}{secrgroup}. Feedback
on the software is also welcome, including suggestions for additional
documentation or new features consistent with the overall design.


`Inverse prediction' uses methods from multivariate calibration (Brown
1982). The goal is to estimate population density (D) and the parameters
of a detection function (usually g0 or lambda0 and sigma) by `matching' statistics
from \code{proxyfn(capthist)} (the target vector) and statistics
from simulations of a 2-D population using the postulated detection
model. Statistics (see Note) are defined by the proxy function,
which should return a vector equal in length to the number of parameters
(default np = 3). Simulations of the 2-D population use either internal C++ code or \code{\link{sim.popn}}.

The simulated population is sampled with internal C++ code, \code{\link{sim.capthist}}, or a user-specified function. Simulations match the detector type (e.g., `single' or `multi') and detector layout specified in traps(capthist), including allowance for varying effort if the layout has a \code{\link{usage}} attribute.

The simulated 2-D distribution of animals is Poisson by default.

Simulations are usually conducted on a factorial experimental design in
parameter space - i.e. at the vertices of a cuboid `box' centred on the
working values of the parameters, plus an optional number of centre
points. 

A multivariate linear model is fitted to predict each vector of simulated
proxies from the known parameter values. The number of simulations at
each design point is increased until a specified precision is reached, 
up to the limit max.nsim.

Once a model with sufficient precision has been obtained, a new working
vector of parameter estimates is `predicted' by inverting the linear
model and applying it to the target vector. A working vector
is accepted as the final estimate when it lies within the box; this
reduces the bias from using a linear approximation to extrapolate a
nonlinear function. If the working vector lies outside the box then a new
design is centred on value for each parameter in the working vector.

Once a final estimate is accepted, further simulations are conducted to
estimate the variance-covariance matrix. These also provide a parametric
bootstrap sample to evaluate possible bias. Set var.nsim = 0
to suppress the variance step.

See Efford et al. (2004) for another description of the method, and
Efford et al. (2005) for an application.

If not provided, the starting values are determined automatically
with the **secr** function \code{\link{makeStart}}.

Linear measurements are assumed to be in metres and density in animals
per hectare (10 000 \eqn{\mbox{m}^2}{m^2}).

If \code{ncores > 1} the \pkg{parallel} package is used to create
processes on multiple cores (see \link{Parallel} for more).

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\references{

Brown, P. J. (1982) Multivariate calibration. \emph{Journal of the Royal
Statistical Society, Series B} \bold{44}, 287--321.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G. (2022) secr: Spatially explicit capture--recapture models. 
  R package version 4.5.4. https://CRAN.R-project.org/package=secr/

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture--recapture: likelihood-based methods. In:
D. L. Thompson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer. Pp. 255--269.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
\emph{Animal Biodiversity and Conservation} \bold{27},
 217--228.

Efford, M. G., Warburton, B., Coleman, M. C. and Barker, R. J. (2005) A
field test of two methods for density estimation. \emph{Wildlife Society
Bulletin} \bold{33}, 731--738.

Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978)
Statistical inference from capture data on closed animal populations.
\emph{Wildlife Monographs} \bold{62}.

}

\keyword{ package }

\seealso{

\code{\link{proxy.ms}}
\code{\link[ipsecr]{ipsecr.fit}},
\code{\link[secr]{secr.fit}},
\code{\link[secr]{capthist}},
\code{\link[secr]{mask}}

}
