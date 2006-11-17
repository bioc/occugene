# cloneFunc.R

# Author:	Oliver Will
# Started:	21 Mar 06
# Modified:	22 Mar 06

################################################################
# Notes:	Functions for transposon clonal libraries. Based
# on cloneFunc.R.
# 
# Dependencies: Needs the function binInsert.
#
# Copyright: Released under the GNU General Public License, v2.
################################################################

checkFormat <- function(anno,clone) {
	# The list of ORFs and clones take different formats.
	# Enforce the correct format.
	# anno <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
	# 			  76,79,85,87,93,94,100,100),byrow=T,ncol=2);
	# clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36);
	if (!(is.matrix(anno))) stop("ORFs not in matrix.");
	if (!(is.vector(clone))) {
		{if ((is.data.frame(clone)&&names(clone)=="position")) {
			clone <- clone$position; 
		}
		else stop("clone is not in correct format.")}
	}
	val <- TRUE;
	return(val);
}

etDelta <- function(d,anno,clone) {
	# Find the Efron and Thisted point estimate.
	# Needs to be in the same folder as binInsertHist.R.
	# d :=	Number of clones to be made.
	# anno :=	ORFs
	# clone :=	Locations of insertions.
	# d <- 10;
	# anno <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
	#			  76,79,85,87,93,94,100,100),byrow=T,ncol=2);
	# clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36);
	checkFormat(anno,clone);
	#source("binInsertHist.R");
	count <- binInsert(clone,anno,TRUE);
	# Get it working for binInsertHist
	# count <- binInsertHist(clone,anno,T);
	count <- c(count,length(clone)-sum(count));
	if (d > length(clone)) stop("d has to be < # of clones");
	# I will have to speed up this function.
	# hist and power functions.
	# ex := expect value
	# v  := variance
	ex <- 0;
	v <- 0;
	for (i in 1:max(count)){
		ex <- (ex+
			   length(count[count==i])*
			   (d/length(clone))^i*
			   (-1)^(i+1));
		v <- (v+
			  length(count[count==i])*
			  (d/length(clone))^(2*i));
	}
	val <- NULL;
	val$expected <- ex;
	val$variance <- v;
	return(val);	
}

fFit <- function(anno,clone,TR=TRUE,b0=0,b1=0,b2=0) {
	# Find the fit.
	# ** Better way to start this function.
	# ** Better way to pass parameters in.
	# Problem with the nonlinear fit function.
	# ** Be able to call generic functions on output.
	#anno <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
	#			  76,79,85,87,93,94,100,100),byrow=T,ncol=2);
	#clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36,
	#		   47,98,58,97,78);
	#TR <- T;
	DEBUG <- TR;
	LARGE <- length(clone)>10000;
	largeMult <- 1000;
	
	# Find the overlap. Eventually give pass the control outside.
	# Worry about last ORF over first later.
	k <- length(anno[,1]);
	overlap <- rep(0,length(anno[,1]));
	indexOverlap <- c((anno[1:(k-1),2]>=anno[2:k,1]),FALSE);
	overlap[indexOverlap] <- 1;
	
	# Put all the switches inside this function **.
	expModel <- function(cumul,TR,b0,b1,b2) {
		# Do the starting values here
		if (DEBUG) {TR <- TRUE;}
		if (all(c(b0==0,b1==0,b2==0))) {
			b0 <- max(cumul$noOrfs);
			b1 <- max(cumul$noOrfs)-1;
			b2 <- 1.0e-1;
			if (LARGE) {b2 <- 1.0e-3;}
			if (DEBUG) {print(paste(b0,b1,b2));}
		}
    	orf.st <- c(b0=b0,b1=b1,b2=b2);
    	orf.fm <- nls(noOrfs ~ b0-b1*exp(-b2*x),data=cumul, 					  start=orf.st,trace=TR);
    	return(orf.fm);
	}
	# Turn m and clone into cumul.
	checkFormat(anno,clone);
	cumExp <- 1:length(clone);
	# Hack to get this to work faster.
	if (LARGE) {cumExp <- 1:trunc(length(clone)/largeMult)};
	tmpHist <- binHist(anno,overlap);
	{if (LARGE) {
		 for (i in 1:trunc(length(clone)/largeMult)) {
		 	cumExp[i] <- binInsertHist(clone[1:(
		 				 i*largeMult)],tmpHist);
		 }
	}
	else {
		for (i in 1:length(clone)) {	
			cumExp[i] <- binInsertHist(clone[1:i],tmpHist);
		}
	}}
	cumul <- t(t(cumExp));
	x <- 1:length(cumul[,1]);
	if (LARGE) {x <- seq(largeMult,trunc(length(clone)/largeMult)*
		 				 largeMult,by=largeMult);}
	if (LARGE) {
		cumul <- c(0,cumul,binInsertHist(clone,tmpHist));
		x <- c(0,x,length(clone));
	}
	if (DEBUG) {print(cumul);}
	cumul <- data.frame(x, noOrfs = cumul);
	if (DEBUG) {
		print(cumul);
		print(NULL);
	}
	val <- expModel(cumul,TR,b0,b1,b2);
	return(val);
}

fCumul <- function(x,b0,b1,b2) {
	# Parametric form of the cumulative occupancy distribution.
	# ** Better way to pas in parameters?
	# b0 <- 3;
	# b1 <- 3;
	# b2 <- 0.01;
	val <- b0-b1*exp(-b2*x);
	return(val);
}

unbiasB0 <- function(anno,
				     clone,
				     iter=1000,
				     seed=NULL,
				     alpha=0.05,
				     TR=TRUE) {
	# Compute the unbiased estimate and confidence interval for 
	# b0 using the parametric bootstrap.
	#anno <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
	#			  76,79,85,87,93,94,100,100),byrow=TRUE,ncol=2);
	#clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36,
	#		   47,98,58,97,78);
	#iter <- 10;
	#seed <- 4;
	#alpha <- 0.05;
	
	# Functions
	simB0 <- function(anno,
					  clone,
					  iter=1000,
					  seed=NULL,
					  TR=FALSE) {
		# Return a list of simulates with no essential genes.
		#anno <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
		#			      76,79,85,87,93,94,100,100),byrow=T,
		#				  ncol=2);
		#clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36,
		#		   47,98,58,97,78);
		#iter <- 100;
		#seed <- 4;
		#TR <- TRUE;
		checkFormat(anno,clone);
		# Monte-carlo.
		if (is.numeric(seed)) set.seed(seed);
		val <- 1:iter;
		n <- length(clone);
		k <- dim(anno)[1];
		for (i in 1:iter) {
			if (TR) print(paste("Iteration: ",i,sep=''));
			tmpClone <- sample(1:max(anno),n,replace=TRUE);
			tmp <- coef(fFit(anno,tmpClone,FALSE))[[1]];
			val[i] <- tmp;
		}
		return(val);
	}
	
	# Body
	checkFormat(anno,clone);
	b0 <- coef(fFit(anno,clone,FALSE))[[1]];
	k <- dim(anno)[1];
	bootB0 <- simB0(anno,clone,iter,seed,TR);
	val <- NULL;
	val$b0 <- k+b0-mean(bootB0);
	sortB0 <- sort(bootB0);
	val$CI <- c(b0+k-sortB0[
				min(ceiling(iter*(1-alpha/2)),iter)],
			    b0+k-sortB0[max(floor(iter*(alpha/2)),1)]);
	# Number of essential genes
	#val$b0 <- k-val$b0;
	#val$CI <- k-val$CI;
	return(val);
}

delta0 <- function(d,anno,clone) {
	# Given the annotation and the list of positions
	# fits b's returns a value of f.
	#m <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
	#			  76,79,85,87,93,94,100,100),byrow=T,ncol=2);
	#clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36,
	#		   47,98,58,97,78);
	n <- length(clone);
	checkFormat(anno,clone);
	tmp <- coef(fFit(anno,clone,FALSE));
	val <- (fCumul(n+d,tmp[1],tmp[2],tmp[3])-
			fCumul(n,tmp[1],tmp[2],tmp[3]));
	return(val[[1]]);
}

unbiasDelta0 <- function(d,
						 anno,
						 clone,
				         iter=1000,
				         seed=NULL,
				         alpha=0.05,
				         TR=TRUE) {
	# Return the unbiased estimator and confidence interval of 
	# delta0.
	#d <- 10;
	#m <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
	#			  76,79,85,87,93,94,100,100),byrow=T,ncol=2);
	#clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36,
	#		   47,98,58,97,78);
	#iter <- 10;
	#seed <- 4;
	#TR <- T;
	
	# Functions
	simDelta0 <- function(d,
						  anno,
						  clone,
						  iter=1000,
						  seed=NULL,
						  TR=FALSE) {
		# Estimate the number of knockouts in the next t ORFs.
		# I am having problems abstracting this for a number 
		# of models.
		#d <- 10;
		#m <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
		#			  76,79,85,87,93,94,100,100),byrow=T,ncol=2);
		#clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36,
		#		   47,98,58,97,78);
		#iter <- 10;
		#seed <- 4;
		#TR <- T;
		checkFormat(anno,clone);
		n <- length(clone);
		# Monte-carlo.
		if (is.numeric(seed)) set.seed(seed);
		val <- 1:iter;
		n <- length(clone);
		k <- dim(anno)[1];
		for (i in 1:iter) {
			if (TR) print(paste("Iteration: ",i,sep=''));
			tmpClone <- sample(1:max(anno),n,replace=TRUE);
			tmp <- delta0(d,anno,tmpClone);
			val[i] <- tmp;
		}
		return(val);
	}
	
	# Body
	checkFormat(anno,clone);
	# Set-up the probability vector
	# Approximate the overlaps as suggested in the chapter
	p <- (anno[,2]-anno[,1]+1)/max(anno);
	p <- p/sum(p);
	p <- c(p,1-sum(p));
	n <- length(clone);
	bootDelta0 <- simDelta0(d,anno,clone,iter,seed,TR);
	val <- NULL;
	val$delta0 <- (eMult(n+d,p)-
				   eMult(n,p)+
				   delta0(d,anno,clone)-
				   mean(bootDelta0));
	sortDelta0 <- sort(bootDelta0);
	val$CI <- c(delta0(d,anno,clone)+
				eMult(n+d,p)-
				eMult(n,p)-
				sortDelta0[min(ceiling(iter*(1-alpha/2)),iter)],
				delta0(d,anno,clone)+
				eMult(n+d,p)-
				eMult(n,p)-
				sortDelta0[max(floor(iter*(alpha/2)),1)]);
	return(val);
}
