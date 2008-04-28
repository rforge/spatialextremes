\name{covariance}
\alias{covariance}
\title{Defines and computes covariance functions}
\description{This function defines and computes several covariance
  function either from a fitted ``max-stable'' model; either by
  specifying directly the covariance parameters.
}
\usage{
covariance(fitted, scale, smooth, cov.mod = "whitmat", dist)
}
\arguments{
  \item{fitted}{An object of class ``maxstab''. Most often this will be
    the output of the \code{\link{fitmaxstab}} function. May be missing
    if \code{scale} and \code{smooth} are given.}
  \item{scale,smooth}{The scale and smooth parameters for  the
    covariance function. May be missing if \code{fitted} is given.}
  \item{cov.mod}{Character string. The name of the covariance
    model. Must be one of "whitmat", "cauchy" or "powexp" for the
    Whittle-Matern, Cauchy and Powered Exponential models.}
  \item{dist}{A numeric vector corresponding to the distance at which
    the covariance function should be evaluated. May be missing.}
}
\value{This function returns the covariance function. Eventually, if
  \code{dist} is given, the covariance function is computed for each
  distance given by \code{dist}.
}
\details{
}
\references{
}
\seealso{
}
\examples{
## 1- Calling covariance using fixed covariance parameters
covariance(scale = 1, smooth = 0.5, cov.mod = "whitmat")
covariance(scale = 1, smooth = 0.5, cov.mod = "whitmat",
  dist = seq(0,5, 0.2))
plot(covariance(scale = 1, smooth = 0.5, cov.mod = "whitmat"),
  from = 0, to = 5, ylab = "Whittle-Matern Covariance Function")

## 2- Calling covariance from a fitted model
require(RandomFields)

##Define the coordinate of each location
n.site <- 30
locations <- matrix(rnorm(2*n.site, sd = sqrt(.2)), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,30, .5), maxstable="extr",
                   n = 30)
ms1 <- t(ms0)

##Now define the spatial model for the GEV parameters
param.loc <- -10 + 2 * locations[,2]
param.scale <- 5 + 2 * locations[,1] + locations[,2]^2
param.shape <- rep(0.2, n.site)
##Transform the unit Frechet margins to GEV 
for (i in 1:n.site)
  ms1[,i] <- param.scale[i] * (ms1[,i]^param.shape[i] - 1) /
  param.shape[i] + param.loc[i]

##Define a model for the GEV margins to be fitted
##shape ~ 1 stands for the GEV shape parameter is constant
##over the region
loc.form <- loc ~ lat
scale.form <- scale ~ lon + I(lat^2)
shape.form <- shape ~ 1

\dontrun{
fitted <- fitmaxstab(ms1, locations, "schlather", loc.form, scale.form,
               shape.form)
covariance(fitted)
covariance(fitted, dist = seq(0,5, 0.2))
}
}
\author{Mathieu Ribatet}
\keyword{htest}