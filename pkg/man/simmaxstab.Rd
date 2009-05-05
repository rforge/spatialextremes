\name{rmaxstab}
\alias{rmaxstab}

\title{Simulation of Max-Stable Random Fields}
\description{This functions generates realisation from a max-stable
  random field.}

\usage{
rmaxstab(n, coord, model = "gauss", cov11 = 1, cov12 = 0, cov22 = 1,
grid = FALSE, ...)
}

\arguments{
  \item{n}{Integer. The number of observations.}
  \item{coord}{A matrix that gives the coordinates of each
  location. Each row corresponds to one location.}
  \item{model}{A character string that gives the max-stable
  model. Currenlty only \code{gauss}, i.e. the Smith model, is
  available.}
  \item{cov11, cov12, cov22}{The parameters for the Smith model.}
  \item{grid}{Logical. Does the coordinates represents grid points?}
  \item{\dots}{Unused.}
}
\value{
  A matrix containing observations from the required max-stable
  model. Each column represents one stations. If \code{grid = TRUE}, the
  function returns an array of dimension nrow(coord) x nrow(coord) x n.
}
\references{
  Schlather, M. (2002) Models for Stationary Max-Stable Random
  Fields. \emph{Extremes} \bold{5}:1,  33--44.
  
  Smith, R. L. (1990) Max-stable processes and spatial
  extremes. Unpublished manuscript.
}
\author{Mathieu Ribatet}
\seealso{\code{\fitmaxstab}}
\examples{
set.seed(8)
x <- seq(0, 10, length = 200)
coord <- cbind(x, x)
data <- rmaxstab(1, coord, cov11 = 9/8, cov22 = 9/8, grid = TRUE)
##We change the margins for visibility
filled.contour(x, x, data[,,1]^0.25)
}
\keyword{distribution}
