\name{Internal}

\alias{proxyfn0}
\alias{proxyfn1}
\alias{proxy.ms}

\alias{simpop}
\alias{simCH}

\title{Internal Functions}

\description{ 

Functions called internally by \pkg{ipsecr}.

}

\usage{

proxy.ms(capthist, model = list(D = ~1, NT = ~1), trapdesigndata = NULL, ...)

proxyfn0(capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...)
proxyfn1(capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...)

simpop(mask, D, N, distribution)
simCH(traps, popn, detectfn, detectpar, NT, noccasions, details = list())
      
}

\arguments{

  \item{capthist}{secr capthist object}
  \item{N.estimator}{character name of closed-population estimator}
  \item{model}{named list of model formulae (see \code{\link{ipsecr.fit}})}
  \item{trapdesigndata}{dataframe with one row for each detector and session}
 
  \item{mask}{ secr mask object}
 \item{D}{numeric density in each mask cell}
 \item{N}{integer number of animals to simulate}
 \item{distribution}{character `even' or `poisson'}
 
 \item{traps}{detector locations as \pkg{secr} traps object}
  \item{popn}{animal locations as \pkg{secr} popn object}
  \item{detectfn}{integer code for detection function (see \link{detectfn})}
  \item{detectpar}{named list of parameter values }
  \item{NT}{numeric hazard of non-target interference at each detector}
  \item{noccasions}{integer number of sampling occasions}
  \item{details}{list with optional additional named arguments}
 
}

\details{

  \code{proxy.ms} is the default proxyfn used by \code{\link{ipsecr.fit}}. 

  \code{simpop} is used by \code{\link{ipsecr.fit}} for popmethod 'internal'. It is faster and simpler than the \pkg{secr} function \code{\link{sim.popn}}.

  \code{simCH} is used by \code{\link{ipsecr.fit}} for CHmethod 'internal'. It is faster and simpler than the \pkg{secr} function \code{\link{sim.capthist}}, and optionally simulates non-target interference.
 
 \code{D} and \code{NT} are matrices with one column per session. 
 
}

\value{

proxy.ms -- a numeric vector of length >= 3 corresponding to proxies for 
    multi-session density and non-target interference models.

proxyfn0, proxyfn1 -- a numeric vector of length 3 corresponding to proxies for 
    population size, capture probability intercept and scale of detection.

simCH -- a single-session capthist object.

}

\seealso{

\code{\link{ipsecr.fit}}

}

\references{

Efford, M. G. (2004) Density estimation in live-trapping studies. \emph{Oikos} \bold{106}, 598--610.

}

\examples{

proxy.ms(captdata)

}