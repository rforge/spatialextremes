\name{condrmaxstab}
\alias{condrmaxstab}
\title{
  Conditional simulation of max-linear random fields
}
\description{
  This function generates (approximate) conditional simulation of unit
  Frechet max-stable random fields from max-linear processes.
}
\usage{
condrmaxstab(n, coord, data.coord, data, cov.mod = "gauss", ..., grid =
FALSE, p = 5000)
}
\arguments{
  \item{n}{Integer. The number of conditional simulations.}
  \item{coord}{A numeric vector or matrix specifying the coordinates
    where the process has to be generated. If \code{coord} is a matrix,
    each row specifies one locations.}
  \item{data.coord}{A numeric vector or matrix specifying the
    coordinates where the process is conditioned.}
  \item{data}{A numeric vector giving the conditioning observations.}
  \item{cov.mod}{A character string specifying the max-stable model. See
    section Details.}
  \item{\dots}{The parameters of the max-stable model. See section
    Details.}
  \item{grid}{Logical. Does \code{coord} defines a grid?}
  \item{p}{An integer. The number of unit Frechet random variables used
    in the max-linear approximation.}
}
\details{
  Any unit Frechet max-stable processes \eqn{\{Z(x)\}}{Z(x)} can be
  approximated by a unit Frechet max-linear process, i.e.,
  \deqn{Z(x) \approx \max_{j=1, \ldots, p} f_j(x) Z_j,}{Z(x) ~ max_{j=1,
      \dots, p} f_j(x) Z_j,}
  where \eqn{f_j}{f_j} are non-negative deterministic functions,
  \eqn{p}{p} is a sufficiently large integer and \eqn{Z_j}{Z_j} are
  independent unit Frechet random variables. Note that to ensure unit
  Frechet margins, the following condition has to be satisfied
  \deqn{\sum_{j=1, \ldots, p} f_j(x) = 1,}{sum_{j=1, \dots, p} f_j(x) =
    1,}
  for all \eqn{x}{x}.

  Currently only the discretized Smith model is implemented for which
  \eqn{f_j(x) = c(p) \varphi(x - u_j ; \Sigma)}{f_j(x) = c(p) phi(x -
  u_j ; Sigma)} where \eqn{\varphi(\cdot; \Sigma)}{phi( . ; Sigma)} is the
  zero mean (multivariate) normal density with covariance matrix
  \eqn{\Sigma}{Sigma}, \eqn{u_j}{u_j} is a sequence of deterministic
  points appropriately chosen and \eqn{c(p)}{c(p)} is a constant
  ensuring unit Frechet margins.
}
\value{
  A matrix containing observations from the required max-stable
  model. Each column represents one stations. If \code{grid = TRUE}, the
  function returns an array of dimension nrow(coord) x nrow(coord) x n.
}
\references{
  Wang, Y. and Stoev, S. A. (2011) Conditional Sampling for Max-Stable
  Random Fields. \emph{Advances in Applied Probability}.
}
\author{
  Mathieu Ribatet
}
\section{Warnings}{
  It may happen that some conditional observations are not honored!
  Sometimes taking a larger \code{p}{p} solves the issue.
}
\seealso{
  \code{\link{rmaxstab}}
}
\examples{
## One dimensional conditional simulations
n.cond.site <- 10
cond.coord <- runif(n.cond.site, -10, 10)
data <- rmaxstab(1, cond.coord, "gauss", var = 3)

x <- seq(-10, 10, length = 250)
cond.sim <- condrmaxstab(5, x, cond.coord, data, var = 3)

matplot(x, t(cond.sim), type = "l", lty = 1, pch = 1)
points(cond.coord, data)

## Two dimensional conditional simulation
cond.coord <- matrix(runif(2 * n.cond.site, -10, 10), ncol = 2)
data <- rmaxstab(1, cond.coord, "gauss", cov11 = 4, cov12 = 0, cov22 = 4)

x <- y <- seq(-10, 10, length = 75)
cond.sim <- condrmaxstab(4, cbind(x, y), cond.coord, data, cov11 = 4,
cov12 = 0, cov22 = 4, grid = TRUE)

par(mfrow = c(2, 2), mar = rep(1, 4))
for (i in 1:4){
image(x, y, log(cond.sim[,,i]), col = heat.colors(64), xaxt = "n", yaxt
= "n", bty = "n")
contour(x, y, log(cond.sim[,,i]), add = TRUE)
text(cond.coord[,1], cond.coord[,2], round(log(data), 2), col = 3)
}
}

\keyword{distribution}
