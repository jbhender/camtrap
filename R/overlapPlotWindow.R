#' Plot overlap or excess within specified window.
#' 
#' Densities are plotted for the entire day, but the overlap or excess are only
#' shown within the requested window.
#'
#' This function is intended to resemble \code{\link{overlapPlot}} from the
#' overlap package with credit to its author, Mike Meredith.
#' However, not all functionality from that function has been 
#' put in place here yet -- rug plotting reamins to do. Currently does not 
#' handle wrapping of window around edges of plot  
#'   
#' @inheritParams overlapEstWindow
#' @param t0,t1 Window times are not scaled internally and should chosen by 
#' the user to reflect \code{xscale}. They are shifter internally to reflect 
#' \code{xcenter}.
#' @param type Character string with \code{overlap} to plot the overlap, 
#' \code{excess} to plot the area between the densities or \code{both} (default)
#' for both. Partial matching okay.
#' @param xscale scale for x-axis: default is 24 for time in hours, 
#' \code{xscale=NA} produces a scale in radians
#' @param xcenter where to center the x-axis: 'noon' (default) or 'midnight'
#' @param extend If not \code{NULL}, the plot is extended by 3 hours on both
#' sides to emphasize the circularity and \code{extend} specifies the color.
#' @param linecol The colors to use for plotting the densities curves.
#' @param overlapcol,excesscol The color(s) to use for plotting the overlap
#' or excess, respectively. 
#' @param mark.window Logical indicating whether to plot vertical lines at the
#' window boundaries.
#' 
#' @return None. Plots to the current device.
#' @export
#' @examples
#' data(simulatedData)
#' overlapPlotWindow(tigerObs,pigObs,t0=18,t1=22)
overlapPlotWindow <-function(A,B,t0=0,t1=24,type=c('both','overlap','excess'),
                             xscale=24,xcenter=c('noon','midnight'),
                             extend=NULL,
                             linetype=c('solid','solid'),
                             linecol=c('black','blue'),
                             linewidth=c(1,1),
                             overlapcol='lightblue',excesscol='lightgrey',
                             mark.window=T,
                             adjust=0.8,kmax=3,n.grid=1e3,...)
{
  type <- match.arg(type,c('both','overlap','excess'))
  isMidnt <- match.arg(xcenter) == 'midnight'
  
	# Estimate and adjust bandwidth #
	bw.A <- getBandWidth(A, kmax = kmax) / adjust
  bw.B <- getBandWidth(B, kmax = kmax) / adjust
  if(is.na(bw.A) || is.na(bw.B))
    stop("Bandwidth estimation failed.")
	
  # Get limits of x-axis based on specified options #
  xsc <- ifelse(is.na(xscale),1,xscale/{2*pi})
  if(is.null(extend)){
    xxRad <- seq(0,2*pi,length=n.grid)
  } else{
    xxRad <- seq(-pi/4,9*pi/4,length=n.grid)
  }
  if(isMidnt){
    xxRad <- xxRad - pi
    if(t0 > 12) t0 <- t0-24
    if(t1 > 12) t1 <- t1-24
  }
  xx <- xxRad*xsc
  # Define density functions #
  dA <- function(x){
    densityFit(A,x,bw.A)	
  }
  dB <- function(x){
    densityFit(B,x,bw.B)
  }
  
  # Define functions of interest #
  fmin <- function(x) pmin(dA(x),dB(x))
  fmax <- function(x) pmax(dA(x),dB(x))
  fdif <- function(x) dA(x) - dB(x)

  densA <- dA(xxRad)/xsc
  densB <- dB(xxRad)/xsc
  
  # Get plotting args #
  dots <- list(...)
  if(length(dots)==1 && class(dots[[1]])=='list'){
    dots <- dots[[1]]
  }
  defaultArgs <- list(main=paste(deparse(substitute(A)), "and",
                                 deparse(substitute(B))), 
                      xlab = 'Time',ylab='Density',xaxt='s',
                      bty='o',type='l',xlim=range(xx),
                      ylim=c(0,max(densA,densB)),
                      las=1
                      )
	useArgs <- modifyList(defaultArgs,dots)
  selPlot <- names(useArgs) %in% c(names(as.list(args(plot.default))),
                                   names(par(no.readonly=TRUE)))
  plotArgs <- useArgs[selPlot]
  plotArgs$x <- 0
  plotArgs$y <- 0
  plotArgs$type <- 'n'
  plotArgs$xaxt <- 'n'
	do.call(plot,plotArgs,quote=TRUE)
	
	# Plot time axis #
  if(useArgs$xaxt!='n'){
    if(!is.na(xscale)){
      if(xcenter=='noon'){
        axis(1,seq(0,24,6),paste(seq(0,24,6),'00',sep=':'))
      } else{
        axis(1,seq(-12,12,6),paste(c(12,18,0,6,12),'00',sep=":"))
      }
    } else{
      if(xcenter=='noon'){
        piScale <- c(expression(0),
                     expression(pi/2),
                     expression(pi),
                     expression(3*pi/4),
                     expression(2*pi))
        axis(1,seq(0,2*pi,pi/2),piScale)
      } else{
        piScale <- c(expression(-pi),
                     expression(-pi/2),
                     expression(0),
                     expression(pi/2),
                     expression(pi))
        axis(1,seq(-pi,pi,pi/2),piScale)
      }
    }
  }  
	# plot the curves #
	if(type=='overlap' || type=='both'){
	  # add a polygon for the overlap in the window only #
	  xx2 <- c(t0,seq(t0,t1,length.out=n.grid),t1)
	  yy <- c(0,fmin(seq(t0,t1,length.out=n.grid)/12*pi),0)/xsc
	  polygon(xx2,yy,col=overlapcol,border=NA)	 
	}
	if(type=='excess' || type=='both'){
	  # add a polygon for the excess
	  xx2 <- c(seq(t0,t1,length.out=n.grid),seq(t1,t0,length.out=n.grid))
	  yy <- c(fmax(seq(t0,t1,length.out=n.grid)/12*pi)/xsc,
	          fmin(seq(t1,t0,length.out=n.grid)/12*pi)/xsc)
	  polygon(xx2,yy,col=excesscol,border=NA)	 	
	}
	if(!is.null(extend)){
	  if(isMidnt){
	    wrap <- c(-pi,pi)*xsc
	  } else{
	    wrap <- c(0,2*pi)*xsc
	  }
	  edge <- par('usr')
	  rect(c(edge[1],wrap[2]),rep(edge[3],2),c(wrap[1],edge[2]),rep(edge[4],2),
	       border=NA,col=extend)
	  box(bty=useArgs$bty)
	}
	lines(xx,densA,lty=linetype[1],col=linecol[1],lwd=linewidth[1])
	lines(xx,densB,lty=linetype[2],col=linecol[2],lwd=linewidth[2]) 
  
	if(mark.window) abline(v=c(t0,t1),lty='dashed',col='darkgrey')
	
  return(invisible(data.frame(x=xx,densityA=densA,densityB=densB)))
}

