\name{fCumul}
\alias{fCumul}
\title{Parametric Function for the Cumulative Occupancy Distribution}
\description{
Returns values for parameterized cumulative occupancy distributions.
}
\usage{
fCumul(x,b0,b1,b2)
}
\arguments{
\item{x}{Point to evaluate}
\item{b0}{Parameter b0}
\item{b1}{Parameter b1}
\item{b2}{Parameter b2}
}
\details{
Function fitted to the cumulative occupancy distribution for a multinomial
distribution. Exponential model := b0-b1*exp(-b2*x).
}
\value{Returns a numeric}
\references{See the book chapter O. Will (**) in **}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
x <- 2
b0 <- 3
b1 <- 3
b2 <- 0.01
val <- fCumul(x,b0,b1,b2)
}
\keyword{models}