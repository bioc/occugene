\name{fFit}
\alias{fFit}
\title{Parametric Fit for the Cumulative Occupancy Distribution}
\description{
Parameterizes the cumulative occupancy distribution.
}
\usage{
fFit(anno,clone,TR=TRUE,b0=0,b1=0,b2=.0)
}
\arguments{
\item{anno}{2-column matrix of annotation}
\item{clone}{vector}
\item{TR}{Report a trace}
\item{b0}{Starting value b0}
\item{b1}{Starting value b1}
\item{b2}{Starting value b2}
}
\details{Fits various parametric functions to the occupancy distribution for 
a multinomial. Using the starting values of b0=0, b1=0, and b2=0 forces the 
function to find starting values for you.}
\value{Returns a object.}
\references{See the book chapter O. Will (**) in **}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
data(sampleAnnotation)
data(sampleInsertions)
anno <- cbind(sampleAnnotation$first,sampleAnnotation$last)
clone <- sampleInsertions$position
TR <- TRUE
fm <- fFit(anno,clone,TR)
}
\keyword{models}