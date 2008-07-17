\name{posterior}
\alias{posterior}

\title{MCMC Sampling of Posterior Distributions}
\description{
  Constructing MCMC samples of prior and posterior distibutions for
  max-stable models.
}
\usage{
posterior(n, init, prior, cov.mod, \dots, psd, burn = 0, thin = 1, adj =
FALSE)
}

\arguments{
  \item{n}{The run-length; the number of sampled vectors
    (excluding \code{init}).}
  \item{init}{Numeric vector giving the initial values for the chain,
  taken to be iteration zero.}
  \item{prior}{An object constructed using \code{\link{prior}}.}
  \item{cov.mod}{A character string specifying the covariance model. See
  \code{\link{fitmaxstab}} for further details. A special case is "none"
  which sample from the prior distribution.}
  \item{\dots}{Arguments to the likelihood. Should include \code{data}
    unless \code{cov.mod} is \code{"none"}.}
  \item{psd}{A vector containing standard deviations for proposal
  distributions.}
  \item{burn}{The burn-in period (an integer); the first \code{burn}
    iterations (including \code{init}) are excluded from the chain.}
  \item{thin}{The thinning interval (an integer); iteration \eqn{k}
    is stored only if \eqn{k} mod \code{thin} is zero (and if \eqn{k}
    greater than or equal to \code{burn}).}
  \item{adj}{Logical. If \code{TRUE}, the posterior distribution is
  "adjusted" to have an appropriate variance.}
}
\value{
  A matrix with \code{1+floor(n/thin)-burn} rows.
  Row labels give the iteration numbers.
  Column labels give parameter names.
  
  An attribute \code{ar} is also returned.
  This is a matrix containing acceptence rates in the first row
  (the number of proposals accepted divided by the number of
  iterations) and ``external rates'' in the second (the number of
  proposals that resulted in a zero likelihood, divided by the
  number of iterations).
}

\seealso{\code{\link{prior}}}
\author{Mathieu Ribatet}
\examples{
## 1- Simulating a Spatial field
require(RandomFields)
n.site <- 20
locations <- matrix(rnorm(2*n.site, sd = sqrt(.2)), ncol = 2)
colnames(locations) <- c("long", "latt")
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,60,.7), maxstable="extr",
                   n = 60)
ms0 <- t(ms0)

## 2- Defining the prior distribution
mean.p <- c(scale = 60, smooth = .7)
sd.p <- diag(c(1000, 1000))
psd <- c(1.73, .007)
prior.p <- prior(mean.p, sd.p)
init <- c(40, 1)

## 3- MCMC stage
mc <- posterior(5000, init, prior.p, "whitmat", data = ms0, fit.marge =
                FALSE, coord = locations, psd = psd, thin = 1, burn = 0)
attributes(mc)$ar

## 4- Using coda for the trace plots
require(coda)
plot(mcmc(mc$mc))
}
\keyword{models}
