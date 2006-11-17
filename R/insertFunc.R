# insertFunc.R

################################################################
#
# R Program
#
# Title:  insertFunc.R
#
# Author: Oliver Will
#
# Date Started:  05 Mar 03
# Modified On:   29 Feb 04
# Last Modified: 22 Mar 06
#
# Purpose: Function to put a list of inserts into the correct
# genes. Returns the cumulative number of hits in the order of
# the vector insert.
#
# Copyright: Released under the GNU General Public License, v2.
#
################################################################        

binHist <- function(orf, 
					overlap=NULL, 
					bp=6264403) {
    # Create end points to use in the hist function. 
    # end.pt := Transition point.
    # orfId  := Number of open reading frame. 0 not.
    # overlapId := 1 if overlap between previous and next. 
    # 	0 if not.
    
    k <- length(orf[,1]);
    # Initialize the vectors.
    if (orf[1,1] == 1) {
        #end.pt <- c(0,orf[1,2]+0.5)
        end.pt <- c(-0.5,orf[1,2]+0.5);
        orfId  <- 1;
        overlapId <- 0;
        if (!is.null(overlap)) {
            if (overlap[1] > 0) {
                #end.pt[2] <- end.pt[2]-overlap[1]
                end.pt[2] <- end.pt[2]-overlap[1]-0.5;
                #end.pt    <- c(end.pt,end.pt[2]+overlap[1])
                end.pt    <- c(end.pt,end.pt[2]+overlap[1]+0.5);
                orfId     <- c(orfId,0);
                overlapId <- c(overlapId,1);
            }
            if (overlap[k] > 0) {
                #end.pt[1] <- end.pt[1]+overlap[k]-1;
                end.pt[1] <- end.pt[1]+overlap[k]-0.5;
            }
        }   
    }
    else {
        end.pt <- 0;

        #end.pt <- c(end.pt,orf[1,1]-1,orf[1,2])
        end.pt <- c(end.pt,orf[1,1]-0.5,orf[1,2]+0.5);
        orfId  <- c(0,1);
        overlapId <- c(0,0);
        if (!is.null(overlap)) {
            if (overlap[1] > 0) {
                #end.pt[3] <- end.pt[3]-overlap[1]
                end.pt[3] <- end.pt[3]-overlap[1]+0.5;
                #end.pt    <- c(end.pt,end.pt[3]+overlap[1]+0.5)
                end.pt    <- c(end.pt,end.pt[3]+overlap[1]+0.5);
                orfId     <- c(orfId,0);
                overlapId <- c(overlapId,1);
            }
            if (overlap[k] > 0) {
                #end.pt[1] <- end.pt[2]+overlap[k]-1
                end.pt[1] <- end.pt[2]+overlap[k]-0.5;
            }
        }  
    }
    # Now do the rest.       
    for (i in 2:k) {
        if (orf[i,1]>end.pt[length(end.pt)]) {
            #end.pt <- c(end.pt,orf[i,1]-1,orf[i,2])
            end.pt <- c(end.pt,orf[i,1]-0.5,orf[i,2])
            orfId  <- c(orfId,0,i)
            overlapId <- c(overlapId,0,0)
            if (!is.null(overlap)) {
                if (overlap[i] > 0) {
                    #end.pt[length(end.pt)] <- 
                    #	end.pt[length(end.pt)]-overlap[i]
                    (end.pt[length(end.pt)] <- 
                    	end.pt[length(end.pt)]-overlap[i]-0.5);
                    #end.pt    <- c(end.pt,
                    #	end.pt[length(end.pt)]+overlap[i])
                    end.pt    <- c(end.pt,
                    	end.pt[length(end.pt)]+overlap[i]+0.5);
                    orfId     <- c(orfId,0);
                    overlapId <- c(overlapId,i);
                }
            }
        }  
        else {
            #end.pt[length(end.pt)-1]
            end.pt[length(end.pt)-0.5];
            #end.pt <- c(end.pt,orf[i,2])
            end.pt <- c(end.pt,orf[i,2]+0.5);
            orfId  <- c(orfId,i);
            overlapId <- c(overlapId,0);
            if (!is.null(overlap)) {
                if (overlap[i] > 0) {
                    #end.pt[length(end.pt)] <- 
                    #	end.pt[length(end.pt)]-overlap[i]
                    (end.pt[length(end.pt)] <- 
                    	end.pt[length(end.pt)]-overlap[i]-0.5);
                    #end.pt    <- c(end.pt,
                    #	end.pt[length(end.pt)]+overlap[i])
                    end.pt    <- c(end.pt,
                    	end.pt[length(end.pt)]+overlap[i]+0.5);
                    orfId     <- c(orfId,0);
                    overlapId <- c(overlapId,i);
                }
            }
        }
    }
    # Take care of the end
    if (orf[k,2]<bp) {
        end.pt <- c(end.pt,bp);
        orfId  <- c(orfId,0);
        overlapId <- c(overlapId,0);
    }
    # Return the ans.
    ans <- NULL;
    ans$end.pt <- end.pt;
    ans$orf    <- orfId;
    ans$overlap <- overlapId;
    return(ans)
}

binInsertHist <- function(insert, 
						  orfHist, 
						  returnCounts=FALSE) {
    # Insert is a list of bp were the inserts occur. 
    # Default is to return the number of different ORFs hit. 
    # orf := End points for a histogram object.
    # returnCounts := Flag for returning the number of times each
    #                 ORF is hit.

    # PARAMETERS
    DEBUG <- FALSE;

    # VARIABLES

    # BODY
    inOrfHist <- 1:length(orfHist$end.pt);
    inOrf     <- inOrfHist[orfHist$orf>0];
    inOverlap <- inOrfHist[orfHist$overlap>0];
    orfOverlap <- unique(orfHist$overlap[orfHist$overlap>0]);
    if (DEBUG) print(c("orfOverlap",length(orfOverlap))); 
    tmp        <- hist(insert,orfHist$end.pt,plot=FALSE);
    if (DEBUG) print(c("Insertions included",sum(tmp$count)));
    if (DEBUG) print(c("Sum count[inOrf]",sum(tmp$count[inOrf])));
    if (DEBUG) print(c("Sum count[inOverlap]",
    	sum(tmp$count[inOverlap])));
    orfCount   <- tmp$count[inOrf];
    overlapCount <- tmp$count[inOverlap];
    orfCount[orfOverlap] <- orfCount[orfOverlap]+overlapCount;
    orfCount[orfOverlap+1] <- orfCount[orfOverlap+1]+overlapCount;
    if (DEBUG) print(c("Sum orfCount",sum(orfCount)));
    if (returnCounts) {
        non.coding <- (length(insert)-sum(tmp$count[inOrf])-
        			   sum(tmp$count[inOverlap]));
        ans <- NULL;
        ans$count <- c(orfCount,non.coding);
        ans$once  <- c(tmp$count[inOrf],non.coding);
        ans$btw   <- rep(0,length(ans$count));
        ans$btw[orfOverlap] <- tmp$count[inOverlap];
    }
    else {
        ans <- length(orfCount[orfCount>0]);
    }
    return(ans);
}

binInsert <- function(insert, 
					  orf, 
					  returnCounts=FALSE, 
					  overlap=NULL,
					  DEBUG=FALSE) {
    # insert is a list of bp were the inserts occur. 
    # orf is a n x 2 matrix of the beginning and ends of the orfs
    # processed so that they are independent of direction.
    # returnCounts is a flag for returning the number of times 
    # each ORF is hit. 
    # Overlap is the number of bp two ORFs overlap.                     
    
    # VARIABLES
    # countOrf :=	No. inserts in each ORF
	# noOrf :=	No. different ORFs hit on insert n
	# countOrfHit := 	No. ORFs hit so far
	# countInsert := 	No. inserts in ORFs
	# countMiss :=	No. inserts that miss ORFs
    countOrf <- rep(0,length(orf[,1])) 
    noOrf <- rep(0,length(insert))    
    countOrfHit <- 0                   
    countInsert <- 0                  
    countMiss   <- 0                  
    
    # BODY
    for (i in 1:length(insert)) {
        if (DEBUG) {
            print(c('i = ',i))
        }
        minIndex <- 1
        maxIndex <- length(orf[,1])
        j <- floor((maxIndex+minIndex)/2)
        last <- 0
        done <- FALSE
        # Binary insert
        while (!done) {
            if (DEBUG) print(c('j = ',j))
            if (insert[i] <= orf[j,2]) {
                if (insert[i] >= orf[j,1]) {
                    done <- TRUE
                    countOrf[j] <- countOrf[j]+1 
                    if (countOrf[j] == 1) {
                        countOrfHit <- countOrfHit+1
                    }
                    # Insert overlap code.
                    if (!is.null(overlap)) {
                        if (insert[i] >= orf[j,2]-overlap[j]) {
                            countOrf[j+1] <- countOrf[j+1]+1
                        }
                        if (countOrf[j+1] == 1) {
                            countOrfHit <- countOrfHit+1
                        }
                    }
                    countInsert <- countInsert+1
                }
                else {
                    if (maxIndex <= minIndex || last == j) {
                        countMiss <- countMiss+1
                        done <- TRUE
                        last <- 0 
                    }
                    else {
                        maxIndex <- j
                        last <- j
                        j <- floor((maxIndex+minIndex)/2)
                    }
                }
            }
            else {
                if (maxIndex <= minIndex || last == j) {
                    countMiss <- countMiss+1
                    done <- TRUE
                    last <- 0
                }
                else {
                    minIndex <- j
                    last <- j
                    j = floor((maxIndex+minIndex)/2)
                }
            }
        }
        noOrf[i] <- countOrfHit
    }
    if (DEBUG) {
    	print(c("Hits = ", countInsert, "Misses = ", countMiss));
    }
    # Figure out how to get countOrf out of this function later.
    # noOrf is the cumulative 
    # Returns the cumulative count of ORFs hit.
    if (!returnCounts){
        temp <- noOrf
    }
    else {
        temp <- countOrf
    }
    return(temp)
}

