\name{plot.maxstab}
\alias{plot.maxstab}
\title{
Model checking of a fitted max-stable model
}
\description{
This function produces several plots to assess the goodness of fit of a
fitted max-stable model.
}
\usage{plot(fitted, sites)}
\arguments{
  \item{fitted}{An object of class \code{maxstab}. Most often, this will
    be the output of \code{\link{fitmaxstab}} or
    \code{\link{lsmaxstab}}.}
  \item{sites}{A vector of integer of length 4 specifying the locations
    to be used for the model checking. If missing, locations will be
    choosen randomly.}
}
\value{
  Several plots.
}
\details{
  The diagonal plots are return level plots. The lower ones are qq-plots
  (on the Gumbel scale) between observed pairwise maxima for each block,
  e.g. year, and the ones obtained by simulations from the fitted
  model. The upper plot compares the fitted extremal coefficient
  functions to semi-empirical estimates from the F-madogram - see
  \code{\link{fmadogram}}. The two remaining plots are the stations
  locations and a qq-plot of blockwise maxima where the block size is 4.
}
\author{Mathieu Ribatet}
\keyword{hplot}
