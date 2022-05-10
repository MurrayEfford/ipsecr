\name{Internal}

\alias{proxyfn0}
\alias{proxyfn1}
\alias{proxyfn2}
\alias{simCH}

\title{Internal Functions}

\description{ 

Functions called internally by \pkg{ipsecr}.

}

\usage{

proxyfn0(capthist, N.estimator =  c("n", "null","zippin","jackknife"))
proxyfn1(capthist, N.estimator =  c("n", "null","zippin","jackknife"))
proxyfn2(capthist, N.estimator =  c("n", "null","zippin","jackknife"))

simCH(traps, popn, detectfn, detectpar, noccasions)

}

\arguments{

  \item{capthist}{secr capthist object}
  \item{N.estimator}{character name of closed-population estimator}
 
  \item{traps}{detector locations as \pkg{secr} traps object}
  \item{popn}{animal locations as \pkg{secr} popn object}
  \item{detectfn}{integer code for detection function (see \link{detectfn})}
  \item{detectpar}{named list of parameter values }
  \item{noccasions}{integer number of sampling occasions}
 
}

\details{

  \code{proxyfn1} is the default proxyfn used by \code{\link{ipsecr.fit}}. 

  \code{proxyfn2} is the default proxyfn used by \code{\link{ipsecr.fit}} when 
  capthist includes nontarget data. 
  
  \code{simCH} is used by \code{\link{ipsecr.fit}} for CHmethod 'internal'. It is faster and simpler than the \pkg{secr} function \code{\link{sim.capthist}}.
  
}

\value{

proxyfn0, proxyfn1 -- a numeric vector of length 3 corresponding to proxies for 
population size, capture probability intercept and scale of detection.

proxyfn2 -- a numeric vector of length 4 corresponding to proxies for 
population size, capture probability intercept, scale of detection and nontarget 
interference.

simCH -- a single-session capthist object.

}

\seealso{

\code{\link{ipsecr.fit}}

}

\references{

Efford, M. G. (2004) Density estimation in live-trapping studies. \emph{Oikos} \bold{106}, 598--610.

}

\examples{

proxyfn1(captdata)

}