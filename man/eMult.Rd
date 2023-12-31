\name{eMult}
\alias{eMult}
\title{Expected Value of the Occupancy Distribution}
\description{
Returns the expected value of the occupancy distribution based on a 
multinomial distribution.
}
\usage{
eMult(n, p, iter=NULL, seed=NULL, experimental=NULL)
}
\arguments{
\item{n}{number of attempts in the multinomial distribution}
\item{p}{probabilities for landing in a specific bin}
\item{iter}{number of iterations used in the Monte-Carlo approximation}
\item{seed}{seed for the random number generator}
\item{experimental}{access to other functions of multinomials}
}
\details{This functions computes the expected value of the occupancy distribution
for a multinomial. In other words, the expected number of bins with 
at least one ball. The experimental argument "oneBall" computes
expected number of bins with exactly one ball and the experimental argument
"nextTo" computes the expected number of bins with one ball next to a bin
with zero balls. Consider any functionality through the experimental 
argument untested.}
\value{Returns a numeric}
\references{
See the book chapter O. Will (**) in ** for specific details about
this package or Johnson, N. L. and Kotz, S. (1977) \emph{Urn Models and Their 
Application: An Approach to Modern Discrete Probability Theory}. 
John Wiley & Sons, New York, NY.
}
\author{Oliver Will \email{owill4@yahoo.com}}
\examples{
n <- 20
p <- c(seq(10,1,-1),47)/100
p <- p/sum(p)
eMult(n,p)
eMult(n,p,iter=1000,seed=4)
}
\keyword{distribution}