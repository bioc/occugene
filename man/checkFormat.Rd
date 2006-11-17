\name{checkFormat}
\alias{checkFormat}
\title{Checks the Format of Annotation and Insertions}
\description{
Checks the format of the annotation and insertions. 
}
\usage{
checkFormat(anno,clone)
}
\arguments{
\item{anno}{2-column matrix of annotation}
\item{clone}{vector}
}
\details{Checks the format of the annotation and insertions list. Annotation 
has to be a matrix of the first and last target in the ORF. Insertions has to 
be a vector. Will stop if not correct format.}
\value{Returns a boolean.}
\references{See the book chapter O. Will (**) in **}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
data(sampleAnnotation)
data(sampleInsertions)
anno <- cbind(sampleAnnotation$first,sampleAnnotation$last)
clone <- sampleInsertions$position
if (checkFormat(anno,clone)) {print("Looks good.");}
}
\keyword{attribute}