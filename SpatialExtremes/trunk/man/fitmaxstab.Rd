\name{fitmaxstab}
\alias{fitmaxstab}
\title{Fits a max-stable process to data}
\description{This function fits max-stable processes to data using
  pairwise likelihood. Two max-stable characterisations are available:
  the Smith and Schlather representations.
}
\usage{fitmaxstab(data,  coord, model = c("smith", "schlather"), cov.mod
= "whitmat", loc.form, scale.form, shape.form, fit.marge = TRUE, \dots,
warn.inf = TRUE, method = "BFGS", std.err.type = "none", corr = FALSE)}
\arguments{
  \item{data}{A matrix representing the data. Each column corresponds to
    one location.}
  \item{coord}{A matrix that gives the coordinates of each
    location. Each row corresponds to one location.}
  \item{model}{A character string that defines the max-stable
    characterisation. Must be one of "smith" or "schlather".}
  \item{cov.mod}{A character string corresponding the the covariance
    model in the Schlather representation. Must be one of "whitmat",
    "cauchy" or "powexp" for the Whittle-Matern, the Cauchy and the
    Powered Exponential covariance family.}
  \item{loc.form, scale.form, shape.form}{R formulas defining the
    spatial linear model for the GEV parameters. May be missing. See
    section details.}
  \item{fit.marge}{Logical. If \code{TRUE} (default), the GEV parameters
    are estimated ``point-wisely'' or using the linear model given by
    \code{loc.form}, \code{scale.form} and \code{shape.form}. If
    \code{FALSE}, observations are supposed to be unit Frechet
    distributed.}
  \item{\dots}{Several arguments to be passed to the \code{\link{optim}}
    function. See section details.}
  \item{warn.inf}{Logical. If \code{TRUE} (default), users are warned if
    the log-likelihood is infinite at starting values.}
  \item{method}{The method used for the numerical optimisation
    procedure. Must be one of \code{BFGS}, \code{Nelder-Mead},
    \code{CG}, \code{L-BFGS-B} or \code{SANN}. See \code{\link{optim}}
    for details.}
  \item{std.err.type}{Character string. Must be one of "none", "score",
    "grad". If "none", standard errors are not computed. Otherwise,
    standard errors are estimated using the sandwich estimates - see
    section Details.}
  \item{corr}{Logical. If \code{TRUE} (non default), the asymptotic
    correlation matrix is computed.}
}
\value{
  This function returns a object of class \code{maxstab}. Such objects
  are list with components:
  \item{fitted.values}{A vector containing the estimated parameters.}
  \item{std.err}{A vector containing the standard errors.}
  \item{fixed}{A vector containing the parameters of the model that
    have been held fixed.}
  \item{param}{A vector containing all parameters (optimised and fixed).}
  \item{deviance}{The (pairwise) deviance at the maximum pairwise
    likelihood estimates.}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}} - for the \code{mle} method.}
  \item{data}{The data analysed.}
  \item{model}{The max-stable characterisation used.}
  \item{fit.marge}{A logical that specifies if the GEV margins were
    estimated.}
  \item{cov.fun}{The estimated covariance function - for the Schlather
    model only.}
  \item{extCoeff}{The estimated extremal coefficient function.}
  \item{cov.mod}{The covariance model for the spatial structure.}
}
\details{
  The function is a wrapper to more specialised functions
  (\code{schlatherfull, schlatherlm, smithfull, smithlm}). As spatial
  data often deal with a large number of location, it is impossible to
  write analytically the joint distribution. Consequently, the fitting
  procedure substitutes the "full likelihood" for the pairwise
  likelihood.

  Let define \eqn{L_{i,j}(x_{i,j}, \theta)}{L_{i,j}(x_{i,j}, theta)} the
  likelihood for site \eqn{i}{i} and \eqn{j}{j}, where \eqn{i = 1,
  \dots, N-1}{i = 1, \dots, N-1}, \eqn{j = i+1, \dots, N}{j=i+1, \dots, 
  N}, \eqn{N}{N} is the number of site within the region and
  \eqn{x_{i,j}}{x_{i,j}} are the joint observations for site \eqn{i}{i}
  and \eqn{j}{j}. Then the pairwise likelihood
  \eqn{PL(\theta)}{PL(theta)} is defined by:

  \deqn{\ell_P = \log PL(\theta) = \sum_{i = 1}^{N-1} \sum_{j=i+1}^{N} \log
  L_{i,j} (x_{i,j}, \theta)}{llik_P = log PL(theta) = sum_{i = 1}^{N-1}
  sum_{j=i+1}^{N} log L_{i,j} (x_{i,j}, theta)}

  As pairwise likelihood is an approximation of the ``full likelihood'',
  standard errors cannot be computed directly by the inverse of the
  Fisher information matrix. Instead, a sandwich estimate must be used
  to account for model mispecification e.g.

  \deqn{\hat{\theta} \sim N(\theta, H^{-1} J H^{-1})}{hat(\theta) ~
    N(\theta, H^{-1} J H^{-1})}
  where \eqn{H} is the Fisher information matrix (computed from the
    misspecified model) and \eqn{J} is the variance of the score
    function.

  \eqn{H} is easily estimated using the observed Hessian matrix given by
  the \code{\link{optim}} function. Estimation of \eqn{J} is much more
  difficult. Currently, we propose two different strategies to estimate
  \eqn{J}.

  \describe{
    \item{\code{grad}}{\eqn{J} is estimated from the gradient e.g. \eqn{J
	= \nabla \ell_P \nabla \ell_P^T}{J = grad(llik_P) grad(llik_P)^T}}
    \item{\code{score}}{\eqn{J} is estimated directly from the variance of
      the observed score function.}
  }
}
\references{
  Cox, D. R. and Reid, N. (2004) A note on pseudo-likelihood constructed
  from marginal densities. \emph{Biometrika} \bold{91}: 729--737.

  Schlather, M. (2002) Models for Stationary Max-Stable Random
  Fields. \emph{Extremes} \bold{5}:1,  33--44.

  Smith, R. L. (1990) Max-stable processes and spatial
  extremes. Unpublished.
}
\examples{
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

##Fit a max-stable process
##  1- using the Schlather representation
\dontrun{
fitmaxstab(ms1, locations, "schlather", loc.form, scale.form,
           shape.form)
}

##  2- using the Smith representation
\dontrun{
fitmaxstab(ms1, locations, "smith", loc.form, scale.form,
           shape.form)}

##  3- without any spatial model for the GEV parameters
## Be careful this could be *REALLY* time consuming
\dontrun{fitmaxstab(ms1, locations, "schlather")}

##  4- Fixing the smooth parameter of the Whittle-Matern family
##     to 0.5 - e.g. considering exponential family
\dontrun{fitmaxstab(ms1, locations, "schlather", smooth = 0.5)}
}
\author{Mathieu Ribatet}
\keyword{htest}