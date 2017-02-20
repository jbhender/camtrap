##' Estimate overlap and excess for a fixed window.
##' 
##' Estimates densities for two sets of sighting times and then
##' returns both the overlap and difference within a specified window.
##' 
##' This functions uses t\code{\link{getBandwidth}} and
##' \code{\link{densityFit}} to estimate density functions for the
##' times provided as \code{A} and \code{B}. The various overlap estimates 
##' described in \code{\link{overlapEst}} are forgone in favor of 
##' adapative quadarature using \code{\link{integrate}}. The arguments
##' \code{t0} and \code{t1} are transformed to radian cscale and then
##' used as the limits of integration, wrapped when appropriate.  
##' 
##' @param A,B Numeric vectors of sighting times in [0,2pi]
##' @param t0,t1 The window over which to estimate the overlap and excess stats
##' @param adjust,kmax See \code{\link{overlapEst}}.
##' 
##' @return A list containing: \describe{
##' \item{overlap}{The area under the minimum of the two densities within
##' the window requested using \code{t0} and \code{t1}.}
##' \item{excess}{Contains the area between the density curves within 
##' the window. If the functions cross within the window, this is the net
##' difference. Differencs are for \eqn{f_A} less \eqn{f_B}.}
##' \item{window}{The requested interval.}
##' }
##' 
##' @author James Henderson based on work by 
##' @export
##' 

overlapEstWindow <- function(A,B,t0=0,t1=24,adjust=0.8,kmax=3){
	
	# rescale window from 24hr scale to [0,2pi)
	t0 <- t0/12*pi
	t1 <- t1/12*pi
		
	# Estimate and adjust bandwidth
	bw.A <- getBandWidth(A, kmax = kmax) / adjust
  bw.B <- getBandWidth(B, kmax = kmax) / adjust

	# Define density functions #
	dA <- function(x){
	   densityFit(A,x,bw.A)	
	}
	dB <- function(x){
	 	densityFit(B,x,bw.B)
	}
	 
	# Define functions of interest
	fmin <- function(x) pmin(dA(x),dB(x))
	fdif <- function(x) dA(x) - dB(x)
	# consider fabs <- function(x) abs(dA(x)-dB(x))
	# if using, should find non-zero regions and integrate
     
	if(t0 < t1){
  	overlap <- integrate(fmin,t0,t1)$value
	  excess <- integrate(fdif,t0,t1)$value	
	} else{
		## wrap estimate overnight
		overlap <- integrate(fmin,t0,2*pi)$value + integrate(fmin,0,t1)$value
		excess  <- integrate(fdif,t0,2*pi)$value + integrate(fdif,0,t1)$value
	}
	
	return(list('overlap'=overlap,'excess'=excess, window=c(t0,t1)))
}

# overlapEstWindow_d <- function(d1,d2,t0=0,t1=24){
	# t0 <- t0/12*pi
	# t1 <- t1/12*pi
	
	# fmin <- function(x) pmin(d1(x),d2(x))
	# fdif <- function(x) d1(x) - d2(x)
	
	# if(t0 < t1){
  	  # overlap <- integrate(fmin,t0,t1)$value
	  # excess <- integrate(fdif,t0,t1)$value	
	# } else{
		# ## wrap estimate overnight
		# overlap <- integrate(fmin,t0,2*pi)$value + integrate(fmin,0,t1)$value
		# excess  <- integrate(fdif,t0,2*pi)$value + integrate(fdif,0,t1)$value
	# }
	
	# out <- list('overlap'=overlap,'excess'=excess,window=c(t0,t1))
	# return(out)
# }
