\name{details}
\alias{details}
\alias{fixedbeta}
\alias{param}

\title{ Detail Specification for ipsecr.fit }

\description{

The function \code{ipsecr.fit} allows many options. Some of these are used
infrequently and have been bundled as a single argument \code{details}
to simplify the documentation. They are described here.
}

\section{Detail components}{

\code{details$contrasts} may be used to specify the coding of factor predictors. The value should be suitable for the 'contrasts.arg' argument of \code{\link[stats]{model.matrix}}. See `Trend across sessions' in \href{https://www.otago.ac.nz/density/pdfs/secr-multisession.pdf}{secr-multisession.pdf} for an example.

\code{details$debug} is used only for debugging. In ordinary use it should not be 
changed from the default (FALSE).

\code{details$distribution} specifies the distribution of the number of
individuals detected \eqn{n}; this may be conditional on the number in the
masked area ("binomial") or unconditional ("poisson").
\code{distribution} affects the sampling variance of the estimated
density. The default is "poisson". The component `distribution' may also
take a numeric value larger than nrow(capthist), rather than "binomial"
or "poisson". 

\code{details$factorial} "full" or "fractional" design. Fractional requires the package **FrF2** (Groemping 2014).

\code{details$FrF2args} list of arguments defining a fractional design

\code{details$ignorenontarget} = TRUE causes the function to ignore
non-target information (the attribute 'nontarget').

\code{details$ignoreusage} = TRUE causes the function to ignore
usage (varying effort) information in the traps component. The default
(\code{details$ignoreusage} = FALSE) is to include usage in the model.

\code{details$popmethod} chooses between internal C++ code, the \pkg{secr} function \code{sim.popn}, 
and a user-provided R function with arguments 'mask', 'D' (density per cell of mask) and 'N' (number of individuals to simulate).

\code{details$CHmethod} chooses between internal C++ code, the \pkg{secr} function \code{sim.capthist}, and a user-provided R function with arguments 'traps', 'popn', 'detectfn', 'detectpar', and 'noccasions'.

Defaults

\itemize{
  \item{boxsize1}{ scalar or vector of length np for size of design}
  \item{boxsize2}{ as for \code{boxsize1}; used from second box onwards}
  \item{centre}{ number of centre points in simulation design}
  \item{min.nsim}{ minimum number of simulations per point }
  \item{max.nsim}{ maximum number of simulations per point }
  \item{dev.max}{ tolerance for precision of points in predictor space }
  \item{var.nsim}{ number of additional simulations to estimate variance-covariance matrix }
  \item{min.nbox}{ minimum number of attempts to `frame' solution }
  \item{max.nbox}{ maximum number of attempts to `frame' solution }
  \item{max.ntries}{ maximum number of attempts at each simulation }
  }
}

\seealso{
  \code{\link{ipsecr.fit} }
}

\references{

Groemping, U. (2014). R Package FrF2 for Creating and Analyzing Fractional Factorial 2-Level Designs. \emph{Journal of Statistical Software}, \bold{56}, 1--56. 
https://www.jstatsoft.org/article/view/v056i01.

}
\keyword{ models }
