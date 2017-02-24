##' Overlap window permutation tests.
##'
##' Assess the signficance of the overlap or excess in a window using a
##' simulation based permutation test.
##' 
##' This function estimates the probability of observing an
##' overlap/excess in a given interval as extreme as the one observed under the
##' null hypothesis that A and B actually come from the same distribution. By
##' default, this hypothesis is tested by simulating new data A' and B' from a
##' density estimated by concatenating A and B.  Setting \code{parametric = F} 
##' will instead give a tradiational permutation test (not yet implemented.)
##' 
##' @inheritParams overlapEstWindow
##' @param nperm The number of simulated samples or permutations to generate.
##' @param two.sided If \code{TRUE} the p-values returned for the excess are 
##' two-sided.P-values for the overlap are always one-sided.
##' @param overlapRef The reference value against which to compare the
##' observed overlap. The default is 1, which makes sense when comparing subsets
##' of the same species.
##' @param excessRef Like \code{overlapRef}, but for the excess.  The default 
##' value is zero, which should make sense in most cases. 
##' 
##' @return An object of class overlapPermObj containing a table of values and
##' information about the function call. 
##' @export
overlapPerm <- function(A,B,t0=0,t1=24,adjust=0.8,nperm=1e3,kmax=3,
                        parametric=T,two.sided=T,overlapRef=1,excessRef=0){
	 obs <- overlapEstWindow(A,B,t0,t1,adjust,kmax)
	 if(parametric){	 	
	    newSamples <- resample(c(A,B),nperm,adjust=adjust,kmax=kmax)
	    Anew <- newSamples[1:length(A),]
	    Bnew <- newSamples[{length(A)+1}:{length(A)+length(B)},]
	    
	    ## Make parallel, check interupt, and/or print progress
      # res is a matrix with row 1 overlap, row 2 excess
	    res <-
	    sapply(1:nperm,function(i){
	    		unlist(overlapEstWindow(Anew[,i],Bnew[,i],t0,t1,adjust,kmax)[1:2])
	    		}
	    	   )
	 } else{
	 	cat('Non-parametric not implemented yet\n')
	 	
	 }

	 ## Compute p-values
	 p <- c()
	 if(two.sided){
	   p[1] <- mean(abs(res[1,]-overlapRef) >= abs(obs$overlap-overlapRef))	
	   p[2] <- mean(abs(res[2,]-excessRef) >= abs(obs$excess-excessRef))
	 } else{
     p[1] <- ifelse({obs$overlap-overlapRef} >= 0,mean(res[1,] >= obs$overlap),
                    mean(res[1,] <= obs$overlap))
     p[2] <- ifelse({obs$excess-excessRef} >= 0, mean(res[2,] >= obs$excess),
                    mean(res[2,] <= obs$excess))
	 }
	 	
	 table <- matrix(c(obs$overlap,overlapRef,p[1],obs$excess,excessRef,p[2]),
	                 2,3,byrow=T)
	 rownames(table) <- c('overlap','excess')
	 colnames(table) <- c('observed','reference','p-value')
	 
	 out <- list(table=table,window=c(t0,t1),nperm=nperm,
	 		parametric=parametric,two.sided=two.sided,overlapRef=overlapRef,
	 		excessRef=excessRef
	 		)
	 class(out) <- c('overlapPermObj','list')
	 return(out)
}