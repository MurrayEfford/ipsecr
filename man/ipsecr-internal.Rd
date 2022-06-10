\name{Internal}

\alias{proxyfn0}
\alias{proxyfn1}
\alias{proxy.nt}
\alias{simCH}

\title{Internal Functions}

\description{ 

Functions called internally by \pkg{ipsecr}.

}

\usage{

proxyfn0(capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...)
proxyfn1(capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...)
proxy.nt(capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...)
proxy.ms(capthist, model, ...)

simCH(traps, popn, detectfn, detectpar, noccasions, details = list())

}

\arguments{

  \item{capthist}{secr capthist object}
  \item{N.estimator}{character name of closed-population estimator}
  \item{model}{named list of model formulae (see \code{\link{ipsecr.fit}})}
 
  \item{traps}{detector locations as \pkg{secr} traps object}
  \item{popn}{animal locations as \pkg{secr} popn object}
  \item{detectfn}{integer code for detection function (see \link{detectfn})}
  \item{detectpar}{named list of parameter values }
  \item{noccasions}{integer number of sampling occasions}
  \item{details}{list with optional additional named arguments}
 
}

\details{

  \code{proxyfn1} is the default proxyfn used by \code{\link{ipsecr.fit}}. 

  \code{proxy.nt} is the default proxyfn used by \code{\link{ipsecr.fit}} when 
  capthist includes nontarget data. 
  
  \code{proxy.ms} can be used to fit multi-session models. 
  
  \code{simCH} is used by \code{\link{ipsecr.fit}} for CHmethod 'internal'. It is faster and simpler than the \pkg{secr} function \code{\link{sim.capthist}}.
  
}

\value{

proxyfn0, proxyfn1 -- a numeric vector of length 3 corresponding to proxies for 
population size, capture probability intercept and scale of detection.

proxy.nt -- a numeric vector of length 4 corresponding to proxies for 
population size, capture probability intercept, scale of detection and nontarget 
interference.

proxy.ms -- a numeric vector of length >= 3 corresponding to proxies for multi-session density models.

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