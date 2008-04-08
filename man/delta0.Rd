\name{delta0}
\alias{delta0}
\title{Number of New Knockouts}
\description{
Point estimate for the number of new ORF knockouts in the next d clones.
}
\usage{
delta0(d,anno,clone)
}
\arguments{
\item{d}{Number of clones to be made}
\item{anno}{2-column matrix of annotation}
\item{clone}{Vector of insertions}
}
\details{
Use the parametric form of the cumulative occupancy distribution to estimate
the number of new ORF knockouts in the next d clones.
}
\value{
A numeric
}
\references{See the book chapter O. Will (**) in **}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
data(sampleAnnotation)
data(sampleInsertions)
anno <- cbind(sampleAnnotation$first,sampleAnnotation$last)
clone <- sampleInsertions$position
delta0(10,anno,clone)
}
\seealso{unbiasDelta0}
\keyword{nonlinear}