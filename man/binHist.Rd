\name{binHist}
\alias{binHist}
\title{Histogram Breakpoints}
\description{
Returns the histogram breakpoints for fast insertion.
}
\usage{
binHist(orf,overlap=NULL,bp=6264403)
}
\arguments{
\item{orf}{2-column matrix of annotation}
\item{overlap}{number position of overlap}
\item{bp}{number of base pairs in genome}
}
\details{
Returns a vector of breakpoints for the binInsertHist function.
}
\value{
\item{end.pt}{Position of last target}
\item{orf}{orfID}
\item{overlap}{Number of targets in overlap}
}
\references{See the book chapter O. Will (**) in **.}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
# **
}
\seealso{binInsertHist}
\keyword{iteration}