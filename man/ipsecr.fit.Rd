\name{ipsecr.fit}
\alias{ipsecr.fit}

\title{

Spatially Explicit Capture--Recapture by Inverse Prediction

}

\description{

Estimate population density by simulation and inverse prediction (Efford 2004; Efford,
Dawson & Robbins 2004). A restricted range of SECR models may be fitted.

}

\usage{

ipsecr.fit(capthist, proxyfn = proxy.ms, model = list(D ~ 1, g0 ~ 1, sigma ~ 1),
    mask = NULL, buffer = 100, detectfn = "HN", binomN = NULL, start = NULL, 
    link = list(), fixed = list(), timecov = NULL, sessioncov = NULL, 
    details = list(), verify = TRUE, verbose = TRUE, ncores = NULL, 
    seed = NULL, ...)

}

\arguments{

  \item{capthist}{secr capthist object including capture data and detector (trap) layout}
  
  \item{proxyfn}{function to compute proxy from capthist for each coefficient (beta parameter)}
  
  \item{model}{list with optional components each symbolically defining a linear 
      predictor for one real parameter using \code{formula} notation }
  
  \item{mask}{\code{\link{mask}} object }
  
  \item{buffer}{scalar mask buffer radius in metres if \code{mask} not specified}

  \item{detectfn}{integer code or character string for shape of
      detection function 0 = halfnormal, 1 = hazard rate etc. -- see \link{detectfn}}
  
  \item{binomN}{integer code for distribution of counts (see Details) }
  
  \item{start}{vector of initial values for beta parameters, or \code{ipsecr} 
      object from which they may be derived }
  
  \item{link}{list with optional components corresponding to `real'
      parameters (e.g., `D', `g0', `sigma'), each a character string in
      \{"log", "logit", "identity", "sin"\} for the link function of one real parameter }
  
  \item{fixed}{list with optional components corresponding to real parameters 
      giving the scalar value to which the parameter is to be fixed }

  \item{timecov}{optional dataframe of values of time (occasion-specific) covariate(s). NOT USED }

  \item{sessioncov}{optional dataframe of values of session-specific covariate(s)}

  \item{details}{list of additional settings, to control estimation (see Details) }

  \item{verify}{logical, if TRUE the input data are checked with \code{\link{verify}} }

  \item{verbose}{ logical, if TRUE then messages are output during execution}
  
  \item{ncores}{ integer number of cores to use for parallel processing}

  \item{seed}{either NULL or an integer that will be used in a call to \code{set.seed}}

  \item{\dots}{ other arguments passed to proxy function }

}

\details{

The vignette should be consulted for a full exposition.

\code{details} is used for various specialized settings listed below. These are 
also described separately - see \code{\link{details}}.

\tabular{lll}{
  Name \tab Default \tab Description \cr
  boxsize1 \tab 0.2  \tab scalar or vector of length np for size of design\cr 
  boxsize2 \tab 0.05  \tab as for \code{boxsize1}; used from second box onwards\cr
  boxtype \tab 'absolute' \tab `absolute' or `relative' \cr
  centre \tab 3 \tab number of centre points in simulation design\cr
  dev.max \tab  0.002 \tab tolerance for precision of points in predictor space \cr
  var.nsim \tab 2000 \tab number of additional simulations to estimate variance-covariance matrix \cr
  keep.sim \tab FALSE \tab if true then the variance simulations are saved \cr
  min.nsim \tab 20 \tab minimum number of simulations per point \cr
  max.nsim \tab 200 \tab maximum number of simulations per point \cr
  min.nbox \tab 2 \tab minimum number of attempts to `frame' solution \cr
  max.nbox \tab 5 \tab maximum number of attempts to `frame' solution \cr
  max.ntries \tab 2 \tab maximum number of attempts at each simulation \cr
  distribution \tab `poisson' \tab `poisson', `binomial' or `even' \cr
  binomN \tab 0 \tab integer code for distribution of counts (unused) \cr
  ignorenontarget \tab FALSE \tab override nontarget attribute of capthist \cr
  ignoreusage \tab FALSE \tab override usage in traps object of capthist \cr
  debug \tab FALSE \tab stop at arbitrary points in execution (varies) \cr
  savecall \tab TRUE \tab optionally suppress saving of call \cr
  newdetector \tab NULL \tab detector type to override detector(traps(capthist)) \cr
  contrasts \tab NULL \tab coding of factor predictors \cr
  popmethod \tab `internal' \tab `internal' or `sim.popn' or a user-provided function \cr
  CHmethod \tab `internal' \tab `internal' or `sim.capthist' or a user-provided function \cr
  factorial \tab `full' \tab `full' or `fractional' design \cr
  FrF2args \tab NULL \tab arguments for FrF2 when factorial = 'fractional' \cr
  extraparam \tab NULL \tab list of starting values for extra parameters (see vignette) \cr

}

}

\value{

An object of class 'ipsecr', a list comprising:

  \item{call}{ the function call (if details$savecall)}
  \item{capthist}{input}
  \item{proxyfn}{input}
  \item{model}{input} 
  \item{mask}{input}
  \item{detectfn}{input}
  \item{start}{input}
  \item{link}{input}
  \item{fixed}{input}
  \item{timecov}{input}
  \item{sessioncov}{input}
  \item{details}{input}
  \item{designD}{list of design data for density}
  \item{trapdesigndata}{list of design data for trap-specific models}
  \item{parindx}{mapping of coefficients (beta parameters) to real parameters}
  \item{vars}{names of covariates in model}
  \item{betanames}{names of coefficients}
  \item{realnames}{names of 'real' parameters}
  \item{code}{integer completion code: 
  1 successful, 
  2 target not within final box, 
  3 exceeded maximum simulations}
  \item{beta}{estimates of coefficients on link scale }
  \item{beta.vcov}{variance-covariance matrix of estimates }
  \item{designbeta}{vertices of final box (design points)}
  \item{sim.lm}{last lm model fit}
  \item{ip.nsim}{total number of simulations }
  \item{var.nsim.OK}{number of successful variance simulations}
  \item{simulations}{ optional simulation output (see details$keep.sim)}
  \item{parameters}{ optional simulation input (see details$keep.sim)}
  \item{variance.bootstrap}{ dataframe summarising simulations for variance estimation }
  \item{version}{package version}
  \item{starttime}{time execution started}
  \item{proctime}{ processor time (seconds) }
  \item{seed}{ RNG state }
  
(The order and composition of the output list may change).

}

\references{

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
\emph{Animal Biodiversity and Conservation} \bold{27},
 217--228.

}

\seealso{

\code{\link{proxy.ms}}, 
\code{\link{predict.ipsecr}}, 
\code{\link{summary.ipsecr}}

}

\examples{

\donttest{

ipsecrdemo <- ipsecr.fit(captdata, ncores = 1, buffer = 100, detectfn = 14, seed = 1237)

}

}
