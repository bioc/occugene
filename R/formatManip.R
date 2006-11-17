# formatManip.R

# Author:	Oliver Will
# Started:	21 Mar 06
# Modified:	27 Mar 06

################################################################
# Notes:	Based readParameters.R. 
#
# Copyright: Released under the GNU General Public License, v2.
################################################################

loadAnnotation <- function(fileName) {
	# Bring the file into the program for use.
	# ** Add the end processing function later.
	# fileName <- "hypothetical_annotation.txt";
	REMOVEENDS <- F;
	FRAC <- 0.2;
	DEBUG <- F; 
    annotation.data <- read.table(fileName,head=T)
    if (REMOVEENDS) {
		stop("Remove ends has not been implemented.");
    }
    if (any(names(annotation.data)!=
    	  c("idNum","first","last","orientation"))) {
    	stop("Format incorrect");
    }
    return(annotation.data);
}

loadInsertions <- function(fileName) {
	# Load the insertion table. 1 column with locations.
	# fileName <- "hypothetical_insertions.txt";
	val <- read.table(fileName,head=T);
	if (names(val)!="position") {stop("Format incorrect.");}
    return(val);
}

occup2Negenes <- function(anno,
			    		  clone,
			    		  INTERGENIC=FALSE) {
	# Convert the occupancy annotation to negenes annotation.
	# Convert the insertion genes list into negenes format.
	# ** The INTERGENIC flag means to encode the intergenic 
	# region
	# I specify in the book chapter.
	#INTERGENIC <- FALSE;
	#anno <- matrix(c(8,17,23,31,37,44,50,56,62,67,66,70,
	# 			    76,79,85,87,93,94,100,100),byrow=TRUE,ncol=2);
	#clone <- c(58,1,25,24,83,23,75,91,96,7,78,25,9,97,36);
	
	#Functions

	# Body
	checkFormat(anno,clone);
	k <- length(anno[,1]);
	n.sites2 <- rep(0,k);
	n.sites <- anno[,2]-anno[,1]+1;
	if (any(n.sites)<=0) {
		stop("Annotation is in the wrong format.");
	}
	n.sites2[1:(k-1)] <- anno[1:(k-1),2]-anno[2:k,1]+1;
	n.sites2[n.sites2<0] <- 0;
	n.sites <- n.sites-n.sites2;
	{if (!INTERGENIC) {
		n.sites[2:k] <- n.sites[2:k]-n.sites2[1:(k-1)];
		n.sites[1] <- n.sites[1]-n.sites2[k];
	}
	else {stop("Intergenic not implemented.")}}
	val <- as.matrix(cbind(n.sites,n.sites2));
	# Convert a list of insertions into Negenes 2 vector format.
	# I need to speed this up and integrate binInsertHist.
	# How fast would this loop be:
	# ** for (i in 1:k) {
		# clone%in%mList[[i]]
	# }
	# Body
	n <- length(clone);
	k <- length(anno[,1]);
	counts <- rep(0,k);
	counts2 <- rep(0,k);
	overlap <- rep(0,k);
	overlap[n.sites2>0] <- 1;
	histTmp <- binHist(anno,overlap);
	countTmp <- binInsertHist(clone,histTmp,TRUE);
	counts <- countTmp$once[1:k];
	counts2 <- countTmp$btw[1:k];
	val <- cbind(val,as.matrix(cbind(counts,counts2)));
	val <- data.frame(val);
	return (val); 
}
