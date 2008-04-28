\name{profile}
\alias{profile}
\alias{profile.maxstab}
\title{Method for profiling fitted max-stable objects}
\description{Computes profile traces for fitted max-stable models.
}
\usage{
\method{profile}{maxstab}(fitted, param, range, n = 10, plot = TRUE,
\dots)
}
\arguments{
  \item{fitted}{An object of class ``maxstab''. Most often, it will be
    the output of the function \code{\link{fitmaxstab}}.}
  \item{param}{A character string giving the model parameter that are to
    be profiled.}
  \item{range}{The range for the profiled model parameter that must be
    explored.}
  \item{n}{Integer. The number of profiled model parameter that must be
    considered.}
  \item{plot}{Logical. If \code{TRUE} (default), the profile trace is
    plotted.}
  \item{\dots}{Extra options that must be passed to the
    \code{\link{plot}} function.}
}
\section{Warnings}{
  This function can be \bold{really} time consuming!
}
\value{
  A matrix. The first column corresponds to the values for which the
  profiled model parameter is fixed. The second column gives the value
  of the pairwise log-likelihood. The remaining columns contain the
  constrained maximum likelihood estimates for the remaining model
  parameters.
}
\examples{
require(RandomFields)

##Define the coordinates of each location
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

##Fit a max-stable process
##  1- using the Schlather representation
\dontrun{
fitted <- fitmaxstab(ms1, locations, "schlather", loc.form, scale.form,
                     shape.form)
}

##Plot the profile pairwise log-likelihood for the smooth parameter
\dontrun{
profile(fitted, "smooth", range = c(.1, 1))
}
}
\author{Mathieu Ribatet}
\keyword{models}