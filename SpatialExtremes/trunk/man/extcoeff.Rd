\name{extcoeff}
\alias{extcoeff}
\title{Plots the extremal coefficient}
\description{
  Plots the extremal coefficient evolution as the coordinates evolves.
}
\usage{
extcoeff(fitted, n = 150, \dots)
}
\arguments{
  \item{fitted}{A object of class "\code{maxstab}". Most often, it will
    be the output of the function \code{\link{fitmaxstab}}}
  \item{n}{Numeric. \code{n^2} corresponds to the total number of
    estimated extremal coefficients for the contour plot.}
  \item{\dots}{Several options to be passed to the \code{\link{contour}}
    function.}
}
\value{
  A plot.
}
\details{
}
\seealso{
  \code{\link{fitmaxstab}}
}
\examples{
## 1- Random field generation
require(RandomFields)
n.site <- 30
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("long", "latt")

ms0 <- MaxStableRF(locations[,1], locations[,2], model = "wh",
                   grid = FALSE, maxstable = "extr",
                   param = c(0, 1, 0, 30, 1), n = 60)
ms0 <- t(ms0)

## 2- Fit 2 max-stable processes
smith <- fitmaxstab(ms0, locations, "gauss",
                    fit.marge = FALSE)
schlather <- fitmaxstab(ms0, locations, "whitmat",
                        fit.marge = FALSE)

## 3- Plot the extremal coefficient
extcoeff(schlather)
extcoeff(smith, add = TRUE, col = 2)
legend("topright", c("Schlather", "Smith"), col = 1:2,
       inset = .05, lty = 1)
}
\author{Mathieu Ribatet}
\keyword{hplot}
