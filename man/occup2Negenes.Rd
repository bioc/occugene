\name{occup2Negenes}
\alias{occup2Negenes}
\title{Convert Occupancy Format to Negenes}
\description{
Convert the annotation and insertion formation of the occupancy package into
the format for the negenes package.
}
\usage{
occup2Negenes(anno,clone,INTERGENIC=FALSE)
}
\arguments{
\item{anno}{2-column matrix of annotation}
\item{clone}{vector of insertion locations}
\item{INTERGENIC}{Process the intergenic region as last ORF.}
}
\details{
Convert the annotation and insertion formation of the occupancy package into
the format for the negenes package. Of the returned data frame, column 1 is 
n.sites, column 2, n.sites2, column 3, counts, column 4, counts2.
}
\value{Returns a data frame}
\references{See the book chapter O. Will (**) in **}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
data(sampleAnnotation)
data(sampleInsertions)
anno <- cbind(sampleAnnotation$first,sampleAnnotation$last)
clone <- sampleInsertions$position
occup2Negenes(anno,clone)
}
\keyword{manip}