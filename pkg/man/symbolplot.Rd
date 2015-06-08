\name{symbolplot}
\alias{symbolplot}
\title{Detecting spatial trends graphically}
\description{This function performs a symbol plot to help in identifying
  any spatial trends}
\usage{
symbolplot(data, coord, which = "gev", plot.border = NULL,
 col = c("#FF000080", "#0000FF80"))
}
\arguments{
  \item{data}{A matrix representing the data. Each column corresponds to
    one location.}
  \item{coord}{A matrix that gives the coordinates of each
    location. Each row corresponds to one location.}
  \item{which}{A character string specifying which values should be
    plotted. Must be one of "gev" (for the GEV marginal parameters),
    "mean" for pointwise mean or "median" for pointwise median.}
  \item{plot.border}{An R function that plots the border of the study
    region. If \code{NULL} (default) no border is plotted.}
  \item{col}{A vector of length 2 giving the colors to be used to fill
    the circles.}
}
\details{
  This function will plot several circles whose center is located at the
  weather stations and whose radius is proportional to the departure of
  the value at that position to the areal mean value.
}
\value{A plot.}
\author{Mathieu Ribatet}
\examples{
## Symbol plot for the Swiss rainfall data set
data(rainfall)
symbolplot(rain, coord, plot.border = swiss)

## Symbol plot the Dutch wind gust data set
data(wind)
data(netherlands)##to produce the border of the Netherlands
symbolplot(wind, coord, plot.border = function() plot(netherlands),
which = "median")
}

\keyword{hplot}
