#' @title  Memory usage of large variable in workspace
#'
#' @description  Memory of variable in decreasing order
#' @param n Number of  variables  needs to be listed ( decreasing memory storage order)
#' @return DataFrame of top memory variables
#' @export
#' @examples
#'
#' largeVariables(n=5)
#'
#'
#'
#'


largeVariables=function(n=3){
  m.df=ll() # return a dataframe that consists of a variable name as rownames, and class and size (in KB) as columns
  #subset(ll(), KB > 1000) # list of object that have over 1000 KB
  Mem.df=m.df[order(m.df$KB),] # sort by the size (ascending)
  Mem.df$GB=Mem.df$KB/10^6
  sum(Mem.df$GB)
  head(Mem.df[order(Mem.df$GB, decreasing=TRUE),],n)

  # oR
  #t=sort( sapply(ls(),function(x){object.size(get(x))}))
  #t/10^9
}
