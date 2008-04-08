\name{unbiasB0}
\alias{unbiasB0}
\title{Unbiased Estimator of the Number of Non-essential ORFs}
\description{
Unbiased point estimate and confidence intervals for the number of non-essential 
ORFs.
}
\usage{
unbiasB0(anno,clone,iter=1000,seed=NULL,alpha=0.05,TR=TRUE)
}
\arguments{
\item{anno}{2-column matrix of annotation}
\item{clone}{Vector of insertions}
\item{iter}{Number of iterations for the bootstrap}
\item{seed}{Seed for the random number generator}
\item{alpha}{Type I error}
\item{TR}{Report a trace}
}
\details{
Fits a parametric function to the cumulative occupancy distribution. Uses
a parametric bootstrap to correct for bias and find confidence intervals for 
the number of non-essential ORFs.
}
\value{
\item{b0}{Unbiased point estimate}
\item{CI}{Confidence interval at the alpha specified}
}
\references{See the book chapter O. Will (**) in **}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
data(sampleAnnotation)
data(sampleInsertions)
anno <- cbind(sampleAnnotation$first,sampleAnnotation$last)
clone <- sampleInsertions$position
TR <- TRUE
iter <- 10
seed <- 4
unbiasB0(anno,clone,iter,seed,TR=TR)
}
\seealso{fFit}
\keyword{nonlinear}