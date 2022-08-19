\name{plot3D.IP}
\alias{plot3D.IP}

\title{

Plot design and saved simulations for one box from a model fitted with \code{\link{ipsecr.fit}}

}

\description{

A 3-D depiction of the design (a box in parameter space) and the resulting simulations (in proxy space).

}

\usage{

plot3D.IP(object, box = 1, oldplot = NULL, plotcentre = TRUE, plotfinal = FALSE, zkludge = -0.2)

}

\arguments{
  \item{object}{ipsecr object from \code{\link{ipsecr.fit}} with details$keep.sim = TRUE}
  \item{box}{integer number of box to plot}
  \item{oldplot}{list containing transofrmations and plot limits from a previous execution}
  \item{plotcentre}{logical; if TRUE the centrepoint of the design box is plotted}
  \item{plotfinal}{logical; if TRUE the final estimates are plotted as a point in parameter space}
  \item{zkludge}{numeric adjustment for base value of z when plotfinal is TRUE}
}

\details{

The function is restricted to single-session models with 3 real parameters. 

A 2-panel plot is generated, so the graphics options should allow at least 2 panels 
(e.g., \code{par(mfrow = c(1,2))}.

Parameters are plotted on the link scale.

The package \pkg{plot3D} is used (Soetaert 2021).

}

\value{

Invisibly returns a list comprising
\item{pmatparm}{pmat used by \pkg{plot3D} for parameter space}
\item{pmatsim}{pmat used by \pkg{plot3D} for proxy space}
\item{pr}{2-row matrix with lower and upper plot limits of each parameter}
\item{sr}{2-row matrix with lower and upper plot limits of each simulated proxy}
}
\seealso{

\code{\link{ipsecr.fit}}

}

\references{

Soetaert, K. (2021). plot3D: Plotting Multi-Dimensional Data. R package version 1.4.
  https://CRAN.R-project.org/package=plot3D
  
}

\examples{
\donttest{
if (requireNamespace("plot3D")) {
    par(mfrow = c(2,2), oma = c(1,1,3,1))
    # plot first box, saving projection and limits for later use
    oldplot <- plot3D.IP(ipsecrdemo, box = 1)
    # plot second box, using projections and limits from first box
    plot3D.IP(ipsecrdemo, box = 2, oldplot, plotfinal = TRUE, zkludge = -0.1)
    mtext(outer = TRUE, side = 3, line = 0.5, adj = c(0.2,0.8), cex = 1.1, 
        c('Parameter space', 'Proxy space'))
}

}

}
