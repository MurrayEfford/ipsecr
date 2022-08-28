\name{details}
\alias{details}
\alias{fixedbeta}
\alias{param}

\title{ Detail Specification for ipsecr.fit }

\description{

The function \code{ipsecr.fit} allows many options. Some of these are used
infrequently and have been bundled as a single argument \code{details}
to simplify the documentation. They are described here in two groups: 
tuning parameters are listed in the following table, and more exotic options 
follow, listed alphabetically.

}

\section{Tuning parameters}{

\tabular{lll}{
Parameter \tab Default \tab Description \cr
boxsize1 \tab 0.2 \tab  scalar or vector of length np for size of design \cr
boxsize2 \tab 0.05 \tab  as for \code{boxsize1}; used from second box onwards \cr
centre \tab 3 \tab    number of centre points in simulation design \cr
min.nsim \tab 20 \tab  minimum number of simulations per point  \cr
max.nsim \tab 200 \tab  maximum number of simulations per point  \cr
dev.max \tab 0.002 \tab   tolerance for precision of points in proxy space  (see below) \cr
var.nsim \tab 2000 \tab  number of additional simulations to estimate variance-covariance matrix  \cr
min.nbox \tab 2 \tab  minimum number of attempts to `frame' solution  \cr
max.nbox \tab 5 \tab  maximum number of attempts to `frame' solution  \cr
max.ntries \tab 2 \tab  maximum number of attempts at each simulation  \cr
}

\code{dev.max} defines a stopping rule: simulations are added in blocks of 
details$min.nsim until a measure of precision is less than dev.max for all proxies. 
If a vector of length 2, the first element applies to the first box and the second to all later boxes. The measure of precision may be the standard error on the link scale (details$boxtype = "absolute") or the coefficient of variation (details$boxtype = "relative").
}

\section{Other `details' components}{

\code{binomN} (default 0 = Poisson) integer code for distribution of counts.

\code{boxtype} (default "absolute") switches between specifying box size as additive on the transformed parameter scale ("absolute") and relative on the transformed parameter scale ("relative").

\code{CHmethod} (default "internal") chooses between internal C++ code, the \pkg{secr} function \code{sim.capthist}, and a user-provided R function with arguments "traps", "popn", "detectfn", "detectpar", and "noccasions". See also popmethod.

\code{contrasts} (default NULL) bmay be used to specify the coding of factor predictors. The value should be suitable for the 'contrasts.arg' argument of \code{\link[stats]{model.matrix}}. See `Trend across sessions' in \href{https://www.otago.ac.nz/density/pdfs/secr-multisession.pdf}{secr-multisession.pdf} for an example.

\code{debug} (default FALSE) is used only for debugging. In ordinary use it should not be 
changed from the default.

\code{distribution} (default "poisson") specifies the distribution of the number of
individuals detected \eqn{n}; this may be conditional on the number in the
masked area ("binomial") or unconditional ("poisson").
\code{distribution} affects the sampling variance of the estimated
density. The component `distribution' may also
take a numeric value larger than nrow(capthist), rather than "binomial"
or "poisson". 

\code{factorial} (default "full") chooses between "full" or "fractional" design. Fractional requires the package **FrF2** (Groemping 2014).

\code{FrF2args} (default NULL) provides a list of arguments defining a fractional design.

\code{ignorenontarget} (default FALSE) may be used to ignore
non-target information (the capthist attribute 'nontarget'). The default is to model non-target information if it is present.

\code{ignoreusage} (default FALSE) may be used to ignore usage (varying effort)
information in the traps component. The default is to include usage in the model.

\code{keep.sim} (default FALSE) when TRUE causes \code{ipsecr.fit} to save additional output, specifically lists (one component per box) of the simulations and parameters for each box, and the final variance simulations.

\code{newdetector} (default NULL) may be used to override the detector type of 
the traps object embedded in the capthist passed to \code{ipsecr.fit}.

\code{Nmax} (default 1e4) maximum allowed population size for simulations.

\code{nontargettype} (default "exclusive") chooses among options "exclusive", "truncated", "erased, "independent", and "dependent". See vignette.

\code{popmethod} (default "internal") chooses between internal C++ code, the \pkg{secr} function \code{sim.popn}, and a user-provided R function with arguments 'mask', 'D' (density per cell of mask) and 'N' (number of individuals to simulate). See also CHmethod.

\code{savecall} (default TRUE) determines whether the function call is included in the output.

\code{YonX} (default TRUE) switches between regression of simulated proxy values Y on controlled parameter values X, or the reverse (which is not fully implemented). 

}

\seealso{
  \code{\link{ipsecr.fit} }
}

\references{

Groemping, U. (2014). R Package FrF2 for Creating and Analyzing Fractional Factorial 2-Level Designs. \emph{Journal of Statistical Software}, \bold{56}, 1--56. 
https://www.jstatsoft.org/article/view/v056i01.

}
\keyword{ models }
