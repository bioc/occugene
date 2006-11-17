\name{binInsert}
\alias{binInsert}
\title{Insert Locations}
\description{
Returns the number of ORF knockouts.
}
\usage{
binInsert(insert,orf,returnCounts=FALSE,overlap=NULL,DEBUG=FALSE)
}
\arguments{
\item{insert}{List of insertion locations}
\item{orf}{2-column matrix of annotation}
\item{returnCounts}{Return the number of insertions}
\item{overlap}{Number of shared targets}
\item{DEBUG}{Flag to debug the code}
}
\details{
Finds the number of ORFs that have an insertion given a list of locations. If 
the returnCounts flag is true, the function returns the number of insertions
per ORF. Uses the function hist for gains in speed.}
\value{Returns a numeric or an object}
\references{See the book chapter O. Will (**) in **.}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
# **
}
\keyword{iteration}