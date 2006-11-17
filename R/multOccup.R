# multOccup.R

################################################################
#
# R Program
#
# Title:  multOccup.R
#
# Author:        Oliver Will
#
# Date Started:  09 Mar 03
# Last Modified: 15 Mar 06
#
# Purpose: Computes the expected value and variance for the
# multinomial occupancy problem.
#
# Dependencies: 
#
# Variables:
# n	:= 	Number of balls
# k	:= 	Number of bins
# p := 	Vector of probabilities for hitting the bins
# iter :=	Number of iterations. If NULL exact computation
#
# Copyright: Released under the GNU General Public License, v2 
# or higher.
#
################################################################

# FUNCTIONS
eMult	<- function(n,p,iter=NULL,seed=NULL,experimental=NULL) {
	# Expected number of occupied bins of a multinomial 
	# distribution. Run by a switch for the internal 
	# functions.
	
	# Internal functions
	equalEMult <- function(n,k) {
		# Expected number of equally probable bins.
		temp <- k*(1-(1-1/k)^n);
		return(temp);
	}	
	equalApproxEMult <- function(n,k,iter=1000,seed=NULL) {
		# Use Monte Carlo itegration.
		if (!(is.null(seed))) {set.seed(seed);}
		tmp <- sample(n,k,replace=TRUE);
		stop("Not implemented yet.");
		return(NULL);
	}
	unequalEMult <- function(n,p) {
		# p vector of probabilities as long as number bins.
		q <- 1-p
		temp <- length(p)-sum(q^n)
		return(temp)
	}
	unequalApproxEMult <- function(n,p,iter=1000,seed=NULL) {
		# p vector of probabilities as long as number of bins.
		# Monte-Carlo integration.
		# I might have to vectorize the loop.
		# n <- 10;
		# p <- (1:5)/sum(1:5);
		# iter <- 100;
		# seed <- NULL;
		if (is.numeric(seed)) set.seed(seed);
		tmp <- matrix(runif(n*iter),nrow=iter);
		tmp2 <- apply(tmp,1,hist,
				   breaks=c(0,cumsum(p)),
				   plot=FALSE);
		occup <- rep(0,iter);
		for (i in 1:iter) {
			tmp3 <- tmp2[[i]];
			occup[i] <-length(tmp3$counts[tmp3$counts>0]);
		}
		return (mean(occup));
	}
	unequalEMult1 <- function(n,p) {
		# Expected number of bins with 1 ball.
		# p vector of probabilities as long as number bins.
		if (sum(p)!=1) {stop("p not a probability vector.");}
		q <- 1-p;
		temp <- n*sum(p*(q)^(n-1));
		return(temp);
	}
	unequalEMult10 <- function(n,p) {
		# Expected number of bins with 1 ball next to a bin 
		# with no balls.
		# p vector of probabilities as long as number bins.
		if (sum(p)!=1) {stop("p not a probability vector.");}
		pPrev <- p[1:(length(p)-1)];
		pNext <- p[2:length(p)];
		qPrev <- 1-pPrev;
		qNext <- 1-pNext;
		v1 <- (qPrev)^(n-1);
		v2 <- pNext*qNext^(n-1);
		temp <- n*(v1%*%v2);
		v1 <- (qNext)^(n-1);
		v2 <- pPrev*qPrev^(n-1);
		temp <- temp+n*(v1%*%v2);
		return(temp);
	}
	
	# Body
	type <- "temp";
	if (is.null(iter)&(length(p)==1)) {type <- "equal";}
	if (!(is.null(iter))&(length(p)==1)) 
		{type <- "approxEqual";}
	if (is.null(iter)&(length(p)!=1)) {type <- "unequal";}
	if (!(is.null(iter))&(length(p)!=1)) 
		{type <- "approxUnequal";}
	if (!(is.null(experimental))) {
		type <- "reset";
		if (experimental=="oneBall") {type <- "oneBall";}
		if (experimental=="nextTo") {type <- "nextTo";}
	}
	val <- switch(type,
		   approxEqual = equalApproxEMult(n,p,iter,seed),
		   equal = equalEMult(n,p),
		   unequal = unequalEMult(n,p),
		   approxUnequal = unequalApproxEMult(n,p,iter,seed),
		   oneBall = unequalEMult1(n,p),
		   nextTo = unequalEMult10(n,p),
		   stop("Parameters not correct")
		   );
	val <- as.numeric(val);
	return(val);
}

varMult	<- function(n,p,iter=NULL,seed=NULL,experimental=NULL) {
	# Variance of occupied bins of a multinomial 
	# distribution. Run by a switch for the internal 
	# functions.
	
	# Functions
	equalVarMult <- function(n,k) {
		temp <- k*(1-1/k)^n
		temp <- temp+k*(k-1)*(1-2/k)^n-k^2*(1-1/k)^(2*n)
		return(temp)
	}
	equalApproxVarMult <- function(n,k) {
		# Use Monte Carlo integration.
		stop("Not implemented yet.");
		return(NULL);	
	}
	unequalVarMult <- function(n,p) {
		# p vector of probabilities as long as number bins.
		temp <- eMult(n,p)-eMult(n,p)^2
		f <- function(i,j) {1-(1-i)^n-(1-j)^n+(1-i-j)^n}
		z <- outer(p,p,f)
		temp <- temp+sum(z)-sum(diag(z))	
		return(temp)
	}
	unequalApproxVarMult <- function(n,p,iter=1000,seed=NULL) {
		# p vector of probabilities as long as number of bins.
		# Monte-Carlo integration.
		# I might have to vectorize the loop.
		# n <- 15;
		# p <- (1:5)/sum(1:5);
		# iter <- 100;
		# seed <- NULL;
		if (is.numeric(seed)) set.seed(seed);
		tmp <- matrix(runif(n*iter),nrow=iter);
		tmp2 <- apply(tmp,1,hist,
				   breaks=c(0,cumsum(p)),
				   plot=FALSE);
		occup <- rep(0,iter);
		for (i in 1:iter) {
			tmp3 <- tmp2[[i]];
			occup[i] <-length(tmp3$counts[tmp3$counts>0]);
		}
		return(var(occup));
	}
	unequalVarMult1 <- function(n,p) {
		# Variance of the number of bins with 1 ball.
		# p vector of probabilities as long as number bins.
		temp <- (eMult(n,p,experimental="oneBall")-
				 eMult(n,p,experimental="oneBall")^2);
		f <- function(i,j) {(n-1)*j*(1-j)^(n-2)*n*i*(1-i)^(n-1)};
		z <- outer(p,p,f);
		temp <- temp+sum(z)-sum(diag(z));	
		return(temp);
	}
	unequalVarMult10 <- function(n,p) {
		# p vector of probabilities as long as number bins.
		# Variance of the number of bins with 1 ball next to a bin 
		# with no balls.
		temp <- (eMult(n,p,experimental="nextTo")-
				 eMult(n,p,experimental="nextTo")^2);
		k <- length(p);
		q <- 1-p;
		prob <- function(n,j) { 
			temp <- (n*p[j+1]*(q[j+1]^(n-1))*(q[j]^n)+n*p[j]*
				(q[j]^(n-1))*(q[j+1]^n));
			return(temp);
		}
		for (i in 1:(k-1)) {
			for (j in 1:(k-1)) {
				# Do the lower and upper triangular part.
				if (abs(i-j) > 1) {
					temp <- temp+prob(n-1,j)*prob(n,i)
				}
				else {
					# Calculate the 3 dependence.
					if (abs(i-j) == 1) {
						l <- min(i,j);
						temp <- (temp+(q[l+2]^(n-1))*
								 n*p[l+1]*(q[l+1]^(n-1))*
								 (q[l]^n));
						temp <- (temp+((n-1)*p[l+2]*(q[l+2]^(n-2))*								 (q[l+1]^(n-1))*
								 (n)*p[l]*(q[l]^(n-1))));
					}
				}
			}
		}
		return(temp);
	}

	# Body
	type <- "temp";
	if (is.null(iter)&(length(p)==1)) {type <- "equal";}
	if (!(is.null(iter))&(length(p)==1)) 
		{type <- "approxEqual";}
	if (is.null(iter)&(length(p)!=1)) {type <- "unequal";}
	if (!(is.null(iter))&(length(p)!=1)) 
		{type <- "approxUnequal";}
	if (!(is.null(experimental))) {
		type <- "reset";
		if (experimental=="oneBall") {type <- "oneBall";}
		if (experimental=="nextTo") {type <- "nextTo";}
	}
	val <- switch(type,
		   approxEqual = equalApproxVarMult(n,p,iter,seed),
		   equal = equalVarMult(n,p),
		   unequal = unequalVarMult(n,p),
		   approxUnequal = unequalApproxVarMult(n,p,iter,seed),
		   oneBall = unequalVarMult1(n,p),
		   nextTo = unequalVarMult10(n,p),
		   stop("Parameters not correct")
		   );
	val <- as.numeric(val);
	return(val);
}