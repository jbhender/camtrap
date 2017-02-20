#' Finds crossings of estimated densities
#' 
#' \code{overlapCross} Fits densities to provided data and returns the 
#' locations where they cross.
#' 
#' This function is designed for use with relatively smooth densities with only
#' a few crossings. It first estimates densities using the provided options
#' and data. It then finds intervals where crossings occur using a grid search 
#' with \code{n.grid} points over [0,2pi] before calling \code{\link{uniroot}}
#' to numerically approximat the crossings in each interval using the
#'difference between densities.
#' 
#' @param A,B Numeric vectors of sighting times in radians, i.e. [0,2pi]
#' @param adjust Scaling factor for the bandwidth estimate, 
#' see \code{\link{overlap::overlapEst}}
#' @param kmax Paramater passed to getBandWidth, see overlapEst
#' @param n.grid The number of grid points used to find rough crossing windows
#' @param scale24 Logical, if true roots are multiplied by 12/pi for 24hr scale 
#' @return A vector of crossings times on the requested scale.
#' 
#' @export
overlapCross <- function(A,B,adjust=0.8,kmax=3,n.grid=128,scale24=T){

  # Return result in radians or hour
  if(scale24){
    scalefac <- 12/pi
  } else{
    scalefac <- 1
  }
  
  grid <- seq(0,2*pi,length=n.grid)
	
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
	f <- function(x){
		dA(x) - dB(x)
	}

	# find approximate roots using grid search
	gA <- dA(grid)
	gB <- dB(grid)
  xx <- which(diff(gA > gB) != 0)
	
  # use uniroot to refine root estimates
  roots <- vector(mode='numeric',length=length(xx))
	for(i in 1:length(xx)){
		roots[i]	 <- uniroot(f,grid[c(0,1)+xx[i]])$root*scalefac
	}

  return(roots)
}
