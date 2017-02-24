##' S3 method print method for overlapPermObj
##' 
##' Tidily print the list created by overlapPerm.
##' 
##' @param obj An object output by \code{\link{overlapPerm}}.
##' @param obj The number of digits to display in the printed table.
##' @export
print.overlapPermObj <- function(obj,digits=3){
  if(!{'overlapPermObj' %in% class(obj)}){
    warning('Object is not of class "overlapPermObj".\n')
  }
  print(round(obj$table,digits))
  sides <- with(obj,ifelse(two.sided,'Two-sided','One-sided'))
  type  <- with(obj,ifelse(parametric,'parametric','non-parametric'))
  
  s <- sprintf('%s permutation tests (%s) for difference from refrence.\n',
               sides,type)
  cat(s)
}