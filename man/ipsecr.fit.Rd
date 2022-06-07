\name{ipsecr.fit}
\alias{ipsecr.fit}

\title{

Spatially Explicit Capture--Recapture by Inverse Prediction

}

\description{

Estimate population density by simulation and inverse prediction (Efford 2004; Efford,
Dawson & Robbins 2004).  A restricted range of SECR models may be fitted.

}

\usage{

ipsecr.fit(capthist, proxyfn = proxyfn1, model = list(D ~ 1, g0 ~ 1, sigma ~ 1),
    mask = NULL, buffer = 100, detectfn = "HN", binomN = NULL, start = NULL, 
    link = list(), fixed = list(), timecov = NULL, sessioncov = NULL, 
    details = list(), verify = TRUE, verbose = TRUE, ncores = NULL, 
    seed = NULL, ...)

}

\arguments{
  \item{capthist}{secr capthist object including capture data and detector (trap) layout}
  \item{proxyfn}{function to compute proxy from capthist for each coefficient (beta parameter)}
  \item{model}{ list with optional components each symbolically defining a linear 
      predictor for one real parameter using \code{formula} notation }
  \item{mask}{ \code{\link{mask}} object }
  \item{buffer}{ scalar mask buffer radius if \code{mask} not specified
      (default 100 m)}

  \item{detectfn}{integer code or character string for shape of
      detection function 0 = halfnormal, 1 = hazard rate etc. -- see \link{detectfn}}
  \item{binomN}{ integer code for distribution of counts (see Details) }
  \item{start}{ vector of initial values for beta parameters, or \code{ipsecr} 
      object from which they may be derived }
  \item{link}{ list with optional components corresponding to `real'
      parameters (e.g., `D', `g0', `sigma'), each a character string in
      \{"log", "logit", "identity", "sin"\} for the link function of one real parameter }
  \item{fixed}{ list with optional components corresponding to real parameters 
      giving the scalar value to which the parameter is to be fixed }

  \item{timecov}{ optional dataframe of values of time (occasion-specific) covariate(s). NOT USED }

  \item{sessioncov}{ optional dataframe of values of session-specific covariate(s). NOT USED }

    \item{details}{ list of additional settings, to control estimation (see Details) }

    \item{verify}{ logical, if TRUE the input data are checked with \code{\link{verify}} }

  \item{verbose}{ logical, if TRUE then messages are output during execution}
    \item{ncores}{ integer number of cores to use for parallel processing}

    \item{seed}{either NULL or an integer that will be used in a call to \code{set.seed}}

  \item{\dots}{ other arguments passed to proxy function }

}

\details{

The vignette should be consulted for a full exposition.

\code{details} is used for various specialized settings listed below. These are
described separately - see \code{\link{details}}.

\tabular{lll}{
  Name \tab Default \tab Description \cr
  boxsize \tab 0.2  \tab scalar or vector of length np for size of design as fraction of central parameter value\cr 
  boxsize2 \tab 0.05  \tab as for \code{boxsize}; used from second box onwards\cr
  boxtype \tab 'absolute' \tab `absolute' or `relative' \cr
  centre \tab 3 \tab number of centre points in simulation design\cr
  dev.max \tab  0.002 \tab tolerance for precision of points in predictor space \cr
  var.nsim \tab 1000 \tab number of additional simulations to estimate variance-covariance matrix \cr
  min.nsim \tab 10 \tab minimum number of simulations per point \cr
  max.nsim \tab 2000 \tab maximum number of simulations per point \cr
  min.nbox \tab 2 \tab minimum number of attempts to `frame' solution \cr
  max.nbox \tab 5 \tab maximum number of attempts to `frame' solution \cr
  max.ntries \tab 2 \tab maximum number of attempts at each simulation \cr
  distribution \tab `poisson' \tab `poisson', `binomial' or `even' \cr
  even \tab FALSE \tab 'binomial' or 'even' population distribution in main simulations\cr
  binomN \tab 0 \tab integer code for distribution of counts (unused) \cr
  ignorenontarget \tab FALSE \tab override nontarget attribute of capthist \cr
  ignoreusage \tab FALSE \tab override usage in traps object of capthist \cr
  debug \tab FALSE \tab stop at arbitrary points in execution (varies) \cr
  savecall \tab TRUE \tab optionally suppress saving of call \cr
  newdetector \tab NULL \tab detector type to override detector(traps(capthist)) \cr
  contrasts \tab NULL \tab coding of factor predictors \cr
  popmethod \tab `internal' \tab `internal' or `sim.popn'\cr
  CHmethod \tab `internal' \tab `internal' or `sim.capthist' or a user-provided function \cr
  factorial \tab `full' \tab `full' or `fractional' design \cr
  FrF2args \tab NULL \tab arguments for FrF2 when factorial = 'fractional' \cr

}

}

\value{

An object of class 'ipsecr', a list comprising:

  \item{call}{ the function call }
  \item{capthist}{input}
  \item{mask}{input}
  \item{detectfn}{input}
  \item{timecov}{input}
  \item{start}{input}
  \item{link}{input}
  \item{fixed}{input}
  \item{model}{input} 
  \item{details}{input}
  \item{designD}{list of design data for density}
  \item{design}{list of design data for detection parameters}
  \item{design0}{list of design data for detection parameters (naive animal)}
  \item{parindx}{mapping of coefficients (beta parameters) to real parameters}
  \item{vars}{names of covariates in model}
  \item{betanames}{names of coefficients}
  \item{realnames}{names of 'real' parameters}
  \item{code}{integer completion code: 1 successful, 2 target not within final box, 3 exceeded maximum simulations}
  \item{beta}{estimates of coefficients on link scale }
  \item{beta.vcov}{variance-covariance matrix of estimates }
  \item{designbeta}{vertices of final box (design points)}
  \item{ip.nsim}{total number of simulations }
  \item{var.nsim.OK}{number of successful variance simulations}
  \item{variance.bootstrap}{ dataframe summarising simulations for variance estimation }
  \item{version}{package version}
  \item{starttime}{time execution started}
  \item{proctime}{ processor time (seconds) }

(The order and composition of the output list is subject to change).

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

\code{\link{proxyfn1}}, 
\code{\link{predict.ipsecr}}, 
\code{\link{summary.ipsecr}}, 
\code{\link{autoini}}

}

\examples{

\donttest{

ipsecrdemo <- ipsecr.fit(captdata, ncores = 1, buffer = 100, detectfn = 14, seed = 1237)

}

}

\section{Warning}{ 

  Simulation becomes unreliable with very sparse
  populations, or sparse sampling, because some simulated datasets will
  have no recaptures or even no captures. The code allows a failed simulation 
  to be repeated (set the `max.ntries' details argument > 1), but results probably 
  should not be relied upon when
  there are warning messages regarding failed simulations. 
  
  }

\note{

Each statistic is expected to have a monotonic relationship with one
parameter when the other parameters are held constant. Typical statistics are -

\tabular{ll}{
Statistic \tab Parameter \cr
\eqn{\hat{N}}{N-hat} \tab \eqn{D} \cr
\eqn{\hat{p}}{p-hat} \tab \eqn{g_0}{g0} \cr
\eqn{RPSV} \tab \eqn{\sigma}{sigma}
}

where \eqn{\hat{N}}{N-hat} and \eqn{\hat{p}}{p-hat} are estimates of population size and capture probability from the
naive application of a nonspatial population estimator, and \eqn{RPSV}
is a trap-revealed measure of the scale of movement.

This method provides nearly unbiased estimates of the detection
parameter g0 when data are from single-catch traps (likelihood-based
estimates of g0 are biased in this case - Efford, Borchers & Byrom 2009).

The implementation largely follows that in Density, and it may help to
consult the Density online help. There are some differences: the M0 and
Mb estimates of population-size in \code{ipsecr.fit} can take non-integer
values; the simulation design used by \code{ipsecr.fit} uses odds(g0)
rather than g0; the default boxsize and CVmax differ from those in
Density 4.5. There is no provision in \code{ipsecr.fit} for two-phase
estimation, using a different experimental design at the second phase.
If you wish you can achieve the same effect by using the estimates as
starting values for a second call of \code{ipsecr.fit} (see examples).

Maximum likelihood estimates from \code{secr.fit} are preferable in
several respects to estimates from inverse prediction (speed*; more complex
models; tools for model selection). \code{ipsecr.fit} is provided for checking
estimates of g0 from single-catch traps, and for historical continuity.

* \code{\link{autoini}} with thin = 1 provides fast estimates
from a simple halfnormal model if variances are not required.

}
