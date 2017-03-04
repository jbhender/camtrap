#' S3 method print method for overlapPermObj
#' 
#' Tidily print the list created by \code{\link{overlapPerm}}.
#' 
#' @param obj An object output by \code{\link{overlapPerm}}.
#' @param obj The number of digits to display in the printed table.
#' @export
print.overlapPermObj <- function(obj,digits=3){
  if(!{'overlapPermObj' %in% class(obj)}){
    warning('Object is not of class "overlapPermObj".\n')
  }
  print(round(obj$table,digits))
  sides <- with(obj,ifelse(two.sided,'two-sided','one-sided'))
  type  <- with(obj,ifelse(parametric,'Parametric','Non-parametric'))
  s <- sprintf('%s permutation tests for difference from refrence.\n',type)
  s1 <- sprintf('For the excess the test is %s.\n',sides)
  cat(s,s1,sep='')
}
