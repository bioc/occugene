\name{unbiasDelta0}
\alias{unbiasDelta0}
\title{Unbiased Number of New Knockouts}
\description{
Unbiased point estimate and confidence intervals for the number of new ORF 
knockouts in the next d clones.
}
\usage{
unbiasDelta0(d,anno,clone,iter=1000,seed=NULL,alpha=0.05,TR=TRUE)
}
\arguments{
\item{d}{Number of new clones}
\item{anno}{2-column matrix of annotation}
\item{clone}{Vector of insertions}
\item{iter}{Number of iterations for the bootstrap}
\item{seed}{Seed for the random number generator}
\item{alpha}{Type I error}
\item{TR}{Report a trace}
}
\details{
Fits a parametric function to the cummulative occupancy distribution. Uses
a parametric bootstrap to correct for bias and find confidence intervals for
the number of new ORF knockouts in the next d clones.
}
\value{
\item{delta0}{Unbiased point estimate}
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
unbiasDelta0(10,anno,clone,iter,seed,TR=TR)
}
\seealso{delta0}
\keyword{nonlinear}