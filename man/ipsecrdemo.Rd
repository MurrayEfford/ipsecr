\name{ipsecrdemo}
\docType{data}

\alias{ipsecrdemo}

\title{ SECR Model Fitted to Demonstration Data }

\description{

  Demonstration data from program Density are provided as a 
  \code{capthist} object (\code{captdata}) ready for input to \code{ipsecr.fit}.

  The fitted models are objects of class \code{ipsecr} formed by 

  \code{ipsecrdemo <- ipsecr.fit(captdata, ncores = 1, detectfn = 14, seed = 1237)}

}

\usage{data(ipsecrdemo)}

\details{

The raw data are 235 fictional captures of 76 animals over 5 occasions
in 100 single-catch traps 30 metres apart on a square grid with origin
at (365,365).

The fitted model uses a hazard halfnormal detection function and default values 
of other arguments.


\tabular{ll}{
  Object \tab Description \cr
  ipsecrdemo \tab fitted ipsecr model -- null \cr
}
 
}

\references{

Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially explicit
capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand. \url{https://www.otago.ac.nz/density/}.

}

\seealso{ \code{\link{capthist}}, \code{\link{read.capthist}}, \code{\link{secrdemo}} }

\examples{

predict(ipsecrdemo)

}

\keyword{datasets}
