\name{gev2frech}
\alias{gev2frech}
\title{Transforms GEV data to unit Frechet ones}
\description{
  Transforms GEV data to unit Frechet ones
}
\usage{
gev2frech(x, loc, scale, shape)
}
\arguments{
  \item{x}{The data to be transformed to unit Frechet}
  \item{loc, scale, shape}{The location, scale and shape parameter of
    the GEV.}
}
\value{
  Returns a numeric vector with the same length of \code{x}
}
\details{
  If Y is a random variable with a GEV distribution with location
  \eqn{\mu}{\mu}, scale \eqn{\sigma}{\sigma} and shape
  \eqn{\xi}{\xi}. Then,

  \deqn{Z = \left[1 + \xi \frac{Y - \mu}{\sigma} \right]_+^{1/\xi}}{Z =
  [1 + \xi (Y - \mu) / \sigma]_+^(1/\xi)}
  is unit Frechet distributed - where \eqn{x_+ = max(0, x)}{x_+ =
    \max(0, x)}.
}
\examples{
x <- c(2.2975896, 1.6448808, 1.3323833, -0.4464904, 2.2737603,
    -0.2581876, 9.5184398, -0.5899699, 0.4974283, -0.8152157)
gev2frech(x, 1, 2, .2)
}
\author{Mathieu Ribatet}
\keyword{datagen}
