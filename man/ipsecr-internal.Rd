\name{Internal}

\alias{proxyfn1}
\alias{proxy.ms}
\alias{detectionDesignData}

\alias{simpop}
\alias{simCH}

\alias{rpsv}
\alias{rpsvi}

\title{Internal Functions}

\description{ 

Functions called internally by \pkg{ipsecr}. 
These are exported and may be called separately for testing.

}

\usage{

proxy.ms(capthist, model = NULL, trapdesigndata = NULL, ...)

detectionDesignData(capthist, byoccasion = FALSE, ...)

proxyfn1(capthist, N.estimator =  c("n", "null","zippin","jackknife"), ...)

simpop(mask, D, N, details = list(), ...)
simCH(traps, popn, detectfn, detparmat, noccasions, NT = NULL, details = 
    list(), ...)
      
rpsv(capthist)
rpsvi(capthist)

}

\arguments{

  \item{capthist}{secr capthist object}
  \item{model}{named list of model formulae (see \code{\link{ipsecr.fit}})}
  \item{trapdesigndata}{dataframe with one row for each detector and session}
  \item{...}{other arguments, mostly unused}
  \item{byoccasion}{logical; if TRUE the output rows are repeated for each occasion }
  
  \item{N.estimator}{character name of closed-population estimator}
 
  \item{mask}{ secr mask object}
  \item{D}{numeric density in each mask cell}
  \item{N}{integer number of animals to simulate}
 
  \item{traps}{detector locations as \pkg{secr} traps object}
  \item{popn}{animal locations as \pkg{secr} popn object}
  \item{detectfn}{integer code for detection function (see \link{detectfn})}
  \item{detparmat}{numeric matrix of detection parameter values}
  \item{noccasions}{integer number of sampling occasions}
  \item{NT}{numeric hazard of non-target interference at each detector}
  \item{details}{list with optional additional named arguments}
 
}

\details{

  \code{proxy.ms} is the default proxyfn used by \code{\link{ipsecr.fit}}. When used internally by \code{\link{ipsecr.fit}}, `model' and `trapdesigndata' are passed automatically. The \dots argument of \code{proxy.ms} may be used to pass arguments to \code{\link{addCovariates}}, especially `spatialdata'. Function \code{detectionDesignData} is used internally to construct design data for non-constant detection models (lambda0, sigma), used in the glm 'data' argument. The capthist argument for \code{detectionDesignData} should always be a list (wrap a single-session capthist in list()).

  \code{simpop} is used by \code{\link{ipsecr.fit}} for popmethod 'internal'. It is faster and simpler than the \pkg{secr} function \code{\link{sim.popn}}. The details component 'distribution' is a character value that may be `poisson' (default) or `even.

  \code{simCH} is used by \code{\link{ipsecr.fit}} for CHmethod 'internal'. It is faster and simpler than the \pkg{secr} function \code{\link{sim.capthist}}, and optionally simulates non-target interference. The argument \code{detparmat} is an individual x parameter matrix, with parameters in the order usual for \code{detectfn}.
 
 \code{D} and \code{NT} are matrices with one column per session. 
 
\code{proxyfn1} is a simple proxy function included mostly for historical reasons. It updates the function of Efford (2004) by log-transforming N, using a complementary log-log transformation instead of odds for p, and using log(RPSV(capthist)) for sigma. If you're interested, look at the code.

\code{rpsv(capthist)} is equivalent to \pkg{secr} RPSV(capthist, CC = TRUE). \code{rpsvi(capthist)} returns a vector of individual-specific rpsv.

}

\value{

proxy.ms -- a numeric vector of length >= 3 corresponding to proxies for a wide range of models, 
    including multi-session density and non-target interference models.
    
detectionDesignData -- a dataframe with one row per individual per session (byoccasion = FALSE) or one row per individual per occasion per session (byoccasion = TRUE), ordered by session, occasion and individual. Columns include x and y coordinates of the individual's centroid, session, and any individual covariates.

proxyfn1 -- a numeric vector of length 3 corresponding to proxies for 
    population size, capture probability intercept and scale of detection.

simpop -- a \code{\link{popn}} object.

simCH -- a single-session \code{\link{capthist}} object.

rpsv -- scalar

rpsvi -- vector, one element per animal

}

\seealso{

\code{\link{ipsecr.fit}}, 
\code{\link{plotProxy}}

}

\references{

Efford, M. G. (2004) Density estimation in live-trapping studies. \emph{Oikos} \bold{106}, 598--610.

}

\note{

\code{proxyfn0} was removed in version 1.2.0.

}

\examples{

proxy.ms(captdata)

}