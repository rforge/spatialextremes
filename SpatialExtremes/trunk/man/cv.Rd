\name{cv}
\alias{cv}
\title{Estimates the penalty coefficient from the cross-validation
  criterion}
\description{
  Estimates the penalty coefficient from the cross-validation
  criterion.
}
\usage{
cv(y, x, knots, degree, pen.range = c(0, 100), \dots)
}
\arguments{
  \item{y}{The response vector.}
  \item{x}{A vector/matrix giving the values of the predictor
    variable(s). If \code{x} is a matrix, each row corresponds to one
    observation.}
  \item{knots}{A vector givint the coordinates of the knots.}
  \item{degree}{The degree of the penalized smoothing spline.}
  \item{pen.range}{A numeric vector giving the lower and upper bound for
    which the cross validation has to be computed.}
  \item{\dots}{Options to be passed to the \code{\link{optimize}}
    function.}
}
\value{
  A list with components 'minimum' and 'objective' which give the
  location of the minimum and the value of the cross-validation
  criterion at that point.
}
\details{
  For every linear smoother e.g. \eqn{\hat{y} = S_\lambda y}{y.hat =
    S_lambda y}, the cross-validation criterion consists in minimizing
  the following quantity:

  \deqn{CV(\lambda) = \sum_{i=1}^n \left(\frac{y_i - \hat{y}_i}{1 -
      S_{\lambda,ii}} \right)^2}{CV(lambda) = \sum_{i=1}^n [(y_i - y.hat_i)
    / (1 - S_{lambda,ii}) ]^2}
  where \eqn{\lambda}{lambda} is the penalty coefficient, \eqn{n}{n} the
      number of observations and \eqn{S_{\lambda,ii}}{S_{lambda,ii}} the
      i-th diagonal element of the matrix \eqn{S_\lambda}{S_lambda}.
}
\references{
  Ruppert, D. Wand, M.P. and Carrol, R.J. (2003) \emph{Semiparametric
    Regression} Cambridge Series in Statistical and Probabilistic
    Mathematics.
}
\seealso{
  \code{\link{cv}}
}
\examples{
n <- 200
x <- runif(n)
fun <- function(x) sin(3 * pi * x)
y <- fun(x) + rnorm(n, 0, sqrt(0.4))
knots <- quantile(x, prob = 1:(n/4) / (n/4 + 1))
cv(y, x, knots = knots, degree = 3)
}
\author{Mathieu Ribatet}
\keyword{htest}