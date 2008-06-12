\name{excofun}
\alias{excofun}
\title{Non parametric estimator of the extremal coefficient function}
\description{Estimates non parametrically the extremal coefficient
  function.
}
\usage{
excofun(data, coord, \dots, prob)
}
\arguments{
  \item{data}{A matrix representing the data. Each column corresponds to
    one location.}
  \item{coord}{A matrix that gives the coordinates of each
    location. Each row corresponds to one location.}
  \item{\dots}{Additional options to be passed to the \code{\link{plot}}
    function.}
  \item{prob}{The probability related to the threshold.}  
}
\details{
  If \code{data} are already componentwise maxima, \code{prob} should be
  zero. Otherwise, users have to define a threshold \eqn{z}{z} (large
  enough) where univariate extreme value arguments are relevant. We
  define \code{prob} such that \eqn{\Pr[Z \leq z] = prob}{Pr[Z <= z] =
    prob}.
}
\value{
  Plots the extremal coefficient function and (invisibly) returns the
  points used for the plot.
}
\references{
  Schlather, M. and Tawn, J. A. (2003) A dependence measure for
  multivariate and spatial extreme values: Properties and
  inference. Biometrika \bold{90}(1):139--156.
}
\seealso{
  \code{\link{madogram}}
}
\examples{
require(RandomFields)
n.site <- 30
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,30, .5), maxstable="extr",
                   n = 40)
ms0 <- t(ms0)
ms1 <- ms0

##Plot the extremal coefficient function
excofun(ms1, locations)
}
\author{Mathieu Ribatet}
\keyword{hplot}
