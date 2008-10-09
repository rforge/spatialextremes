\name{prior}
\alias{prior}

\title{Construction of Prior Distributions}
\description{
  Constructing prior distibutions for a max-stable model.
  parameter.
  
}
\usage{
prior(mean, cov)
}

\arguments{
  \item{mean, cov}{The prior distibution for the parameters of the
  max-stable parameteres multinormal, with mean \code{mean} (a numeric
  vector) and covariance matrix \code{cov} (a symmetric positive
  definite n by n matrix).} 
}
\value{
  Returns an object which is essentially just a list of the arguments
  passed.
}

\seealso{}

\examples{
mat <- diag(c(10000, 10000, 100))
prior(mean = c(0,0,0), cov = mat)
}
\keyword{models}
