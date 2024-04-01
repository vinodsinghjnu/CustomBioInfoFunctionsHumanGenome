#' @title  Percent overlap of a GenomicRange with other
#'
#' @description Percent overlap of a genomic range with other
#' @param FirstContext GenomicRange object. (query: overlapped)
#' @param SecondContext GenomicRange object. (overlapped to)
#' @return Percent overlap of first genomic range to second genomic range
#' @export
#' @examples
#' gr1 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),  ranges = IRanges(start = c(1,110,105), end = c(100, 120, 150 )))
#' gr2 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+", "+"),  ranges = IRanges(start = c(1,12,60, 105), end = c(25, 50, 70, 115 )))
#' pctOverlap_Of_FirstGrToSecondGr(FirstContext=gr1, SecondContext=gr2)
#'
#'
#'
pctOverlap_Of_FirstGrToSecondGr=function(FirstContext, SecondContext){

  (sum(width(GenomicRanges::intersect(FirstContext,SecondContext, ignore.strand=TRUE)))/sum(width(GenomicRanges::reduce(FirstContext))))*100
}
