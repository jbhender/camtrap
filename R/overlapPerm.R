##' Overlap window permutation tests.
##'
##' Assess the signficance of the overlap or excess in a window using a
##' simulation based permutation test.
##' 
##' Specifically, this function estimates the probability of observing an
##' overlap/excess in a given interval as extreme as the one observed under the
##' null hypothesis that A and B actuall come from the same distribution. By
##' default, this hypothesis is tested by simulating new data A' and B' from a
##' density estimated by concatenating A and B.  Setting \code{parametric = F} 
##' will instead give a tradiational permutation test (not yet implemented.)
##' 
##' @inheritParams overalpEstWindow
##' @param nperm The number of simulated samples or permutations to generate.
##' @param two.sided If \code{TRUE} the p-values returned for the excess are 
##' two-sided.P-values for the overlap are always one-sided.

overlapPerm <- function(A,B,t0=0,t1=24,adjust=0.8,nperm=1e3,kmax=3,
                        parametric=T,two.sided=T){
	 obs <- overlapEstWindow(A,B,t0,t1,adjust,kmax)
	 if(parametric){	 	
	    newSamples <- resample(c(A,B),nperm,adjust=adjust,kmax=kmax)
	    Anew <- newSamples[1:length(A),]
	    Bnew <- newSamples[{length(A)+1}:{length(A)+length(B)},]
	    
	    ## Make parallel, check interupt, and/or print progress
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
	   p[1] <- mean(pmin(res[1,],1-res[1,]) <= obs$overlap)	
	   p[2] <- mean(abs(res[2,]) >= abs(obs$excess))
	 } else{
	 	if(obs$overlap > .5){
	 		p[1] <- mean(res[1,] >= obs$overlap)
	 	} else{
	 		p[1] <- mean(res[1,] <= obs$overlap)
	 	}
	 	
	 	if(obs$excess > 0){
	 		p[2] <- mean(res[2,] >= obs$excess)
	 	} else{
	 		p[2] <- mean(res[2,] <= obs$excess)

	 	}
	 }
	 
	 table <- matrix(c(obs$overlap,p[1],obs$excess,p[2]),2,2,byrow=T)
	 rownames(table) <- c('overlap','excess')
	 colnames(table) <- c('observed','p-value')
	 
	 out <- list(table=table,window=c(t0,t1),nperm=nperm,
	 		parametric=parametric,two.sided=two.sided
	 		)
	 class(out) <- c('overlapPermObj')
	 return(out)
}

print.overlapPermObj <- function(obj){
	print(obj$table)
	sides <- ifelse(two.sided,'two-sided','one-sided')
	type  <- ifelse(parametric,'parametric','non-parametric')
	s <- sprintf('')
}
