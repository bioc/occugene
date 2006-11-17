\name{binInsertHist}
\alias{binInsertHist}
\title{Insert Locations Quickly}
\description{
Given a list of locations, returns the number of ORFs hit.
}
\usage{
binInsertHist(insert,orfHist,returnCounts=FALSE)
}
\arguments{
\item{insert}{List of insertion locations}
\item{orfHist}{Histogram breakpoints}
\item{returnCounts}{Return the number of insertions}
}
\details{
Finds the number of ORFs that have an insertion given a list of locations. If 
the returnCounts flag is true, the function returns the number of insertions
per ORF. Uses the function hist for gains in speed.
}
\value{Returns a numeric or an object}
\references{See the book chapter O. Will (**) in ** }
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
# **
}
\seealso{binHist}
\keyword{iteration}