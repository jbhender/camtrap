#' Plot overlap or excess within specified window.
#' 
#' Densities are plotted for the entire day, but the overlap or excess are only
#' shown within the requested window.
#' 
#' @inheritParams overlapEstWindow
#' @param type Character string with \code{overlap} to plot the overlap, 
#' \code{excess} to plot the area between the densities. Partial matching okay.
#' @param mark.window Logical indicating whether to plot vertical lines at the
#' window boundaries.
#' 
#' @return None. Plots to the current device.
#' @export
overlapPlotWindow <-function(A,B,t0=0,t1=24,adjust=0.8,kmax=3,
                             type=c('overlap','excess'),mark.window=T,
                             n.grid=1e3,...)
{
  type <- match.arg(type,c('overlap','excess'))
  
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
	fmax <- function(x) pmax(dA(x),dB(x))
	fdif <- function(x) dA(x) - dB(x)
	 
	# plot the curves #
  curve(dA,0,2*pi,las=1,ylab='density',xlab='time',...)	 
  curve(dB,0,2*pi,add=T,col='blue')	
	 
	if(type=='overlap'){
	  # add a polygon for the overlap in the window only #
		xx <- c(t0,seq(t0,t1,length.out=n.grid),t1)
		yy <- c(0,fmin(seq(t0,t1,length.out=n.grid)),0)
		polygon(xx,yy,col='grey',border=NA)	 
	}
	if(type=='excess'){
		# add a polygon for the excess
		xx <- c(seq(t0,t1,length.out=n.grid),seq(t1,t0,length.out=n.grid))
		yy <- c(fmax(seq(t0,t1,length.out=n.grid)),fmin(seq(t1,t0,length.out=n.grid)))
    polygon(xx,yy,col='lightblue',border=NA)	 	
	}
	if(mark.window) abline(v=c(t0,t1),lty='dashed',col='darkgrey')
}
