\name{etDelta}
\alias{etDelta}
\title{Number of New ORF Knockouts}
\description{
Estimates the number of new knockouts in next d clones.
}
\usage{
etDelta(d,anno,clone)
}
\arguments{
\item{d}{number of new clones}
\item{anno}{2-column matrix of annotation}
\item{clone}{vector}
}
\details{
Estimates the number of new ORF knockouts in the next d clones using the 
method outlined by Efron and Thisted.
}
\value{
\item{expected}{Expected value}
\item{variance}{Variance}
}
\references{See the book chapter O. Will (**) in ** and also Efron, B. and 
Thisted, R. (1976) Estimating the number of unseen species: How many words did 
Shakespere know? \emph{Biometrika}. 63, 435-447.}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
data(sampleAnnotation)
data(sampleInsertions)
a.data <- sampleAnnotation
experiment <- sampleInsertions
orf <- cbind(a.data$first,a.data$last)
clone <- experiment$position
etDelta(10,orf,clone)
}
\keyword{univar}