\name{Internal}

\alias{proxyfn0}
\alias{proxyfn1}

\title{Internal Functions}

\description{ 

Functions called internally by \pkg{ipsecr}.

}

\usage{

proxyfn0(capthist, N.estimator =  c("n", "null","zippin","jackknife"))
proxyfn1(capthist, N.estimator =  c("n", "null","zippin","jackknife"))

}

\arguments{

  \item{capthist}{secr capthist object}
  \item{N.estimator}{character name of closed-population estimator}
  
}

\details{

  \code{proxyfn1} is the default proxyfn used by \code{\link{ipsecr.fit}}. 
  
}

\value{

proxy -- a numeric vector of length 3 corresponding to proxies for population size, 
capture probability intercept and scale of detection.

}

\seealso{

\code{\link{ipsecr.fit}}

}

\references{

Efford, M. G. (2004) Density estimation in live-trapping studies. \emph{Oikos} \bold{106}, 598--610.

}

\examples{

proxyfn0(captdata)

}