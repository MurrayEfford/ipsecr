\name{predict.ipsecr}
\alias{predict.ipsecr}

\title{ SECR Model Predictions }

\description{
Evaluate a spatially explicit capture--recapture model. That is, compute the `real' parameters corresponding to the `beta' parameters of a fitted model for arbitrary levels of any variables in the linear predictor.
}

\usage{

\method{predict}{ipsecr}(object, newdata = NULL, type = c("response", "link"), 
    se.fit = TRUE, alpha = 0.05, savenew = FALSE, \dots)

}

\arguments{
  \item{object}{ \code{ipsecr} object output from \code{ipsecr.fit} }
  \item{newdata}{ optional dataframe of values at which to evaluate model }
  \item{type}{ character; type of prediction required. The default ("response") provides estimates of the `real' parameters.} 
  \item{se.fit}{ logical for whether output should include SE and confidence intervals }
  \item{alpha}{ alpha level for confidence intervals }
  \item{savenew}{ logical for whether newdata should be saved }
  \item{\dots}{ other arguments passed to \code{\link{newdata}} }
}
\details{
  
The variables in the various linear predictors are described in
\href{https://www.otago.ac.nz/density/pdfs/secr-models.pdf}{secr-models.pdf} and listed for the particular model in the
\code{vars} component of \code{object}.

Optional \code{newdata} should be a dataframe with a column for each of
the variables in the model (see `vars' component of \code{object}).  If
\code{newdata} is missing then a dataframe is constructed automatically.

Default \code{newdata} are for a naive animal on the first occasion;
numeric covariates are set to zero and factor covariates to their base
(first) level. The argument `all.levels' may be passed 
to \code{newdata}; if TRUE then the default newdata includes 
all factor levels.

\code{realnames} may be used to select a subset of parameters.

Standard errors for parameters on the response (real) scale are by the
delta method (Lebreton et al. 1992), and confidence intervals are
backtransformed from the link scale.

The value of \code{newdata} is optionally saved as an attribute.

}

\value{
  When \code{se.fit} = FALSE, a dataframe identical to \code{newdata} except for the addition of one column for each `real' parameter. Otherwise, a list with one component for each row in \code{newdata}. Each component is a dataframe with one row for each `real' parameter (density, g0, sigma, b) and columns as below

\tabular{ll}{
link        \tab link function \cr 
estimate    \tab estimate of real parameter \cr
SE.estimate \tab standard error of the estimate \cr
lcl         \tab lower 100(1--alpha)\% confidence limit \cr
ucl         \tab upper 100(1--alpha)\% confidence limit \cr
}

When \code{newdata} has only one row, the structure of the list is
`dissolved' and the return value is one data frame.

For \code{detectpar}, a list with the estimated values of detection
parameters (e.g., g0 and sigma if detectfn = "halfnormal"). In the case
of multi-session data the result is a list of lists (one list per
session).

}

\note{

  \code{\link{predictDsurface}} should be used for predicting density at many
  points from a model with spatial variation. This deals automatically
  with scaling of x- and y-coordinates, and is much faster than
  predict.ipsecr. The resulting Dsurface object has its own plot method.
  
  The argument `scaled' was removed from both predict methods in version 2.10 as the scaleg0 and scalesigma features had been superceded by other parameterisations.
  
}

\references{ 
Lebreton, J.-D., Burnham, K. P., Clobert, J. and Anderson, D. R. (1992) Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. \emph{Ecological Monographs} \bold{62}, 67--118. }

\seealso{ \code{\link{ipsecr.fit}}, \code{\link{predictDsurface}} }

\examples{

predict (ipsecrdemo)

}

\keyword{ models }

