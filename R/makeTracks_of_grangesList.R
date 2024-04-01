#' @title  make tracks list from GenomicRange lists.
#'
#' @description Make tracks from genomic ranges objects list.
#' @param grlist list of genomic range objects.
#' @param location a list of genomic location with fields ´chr´, ´from´ and ´to´. Optional: Must be given if you choose ´if_plain=FALSE´
#' @param assmblyName human genome assembly name (hg19 or hg38 or t2t). Default: hg19
#' @param if_plain : ´TRUE´ then it will return only tracks of given genomic range. Default: FALSE: Tracks of the given chromosomes will also be returned.
#' @return list of tracks for input list of genomic ranges objects.
#' @export
#' @examples
#'
#'
#' data(cpgIslands)
#' mygrlist=list(gr1=cpgIslands,gr2=cpgIslands)
#' loc=list(chr='chr7',from=26700000, to=26750000)
#' tracks=makeTracks_of_grangesList(grlist=mygrlist, location=loc, assmblyName='hg19')
#' plotTracks(tracks, from = loc$from, to = loc$to)
#'
#' # Example 2, only genomic ranges no annotation
#'
#' gr1 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),  ranges = IRanges(start = c(1,110,105), end = c(100, 120, 150 )))
#' gr2 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+", "+"),  ranges = IRanges(start = c(1,12,60, 105), end = c(25, 50, 70, 115 )))
#' plotTracks(makeTracks_of_grangesList(list(gr1=gr1,gr2=gr2), if_plain=TRUE), shape="box")
#'
makeTracks_of_grangesList=function(grlist, location, assmblyName='hg19',if_plain=FALSE){
  if(if_plain==TRUE & missing(location)){
    grTracks.List=purrr::map(grlist, AnnotationTrack )
    for(x in seq_along(names(grlist))) {grTracks.List[[x]]@name = names(grlist)[x]}
    return(grTracks.List)
  }

  if(assmblyName=='hg19'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    BS_genome=BSgenome.Hsapiens.UCSC.hg19
    gen=unique(genome(BS_genome))
  }else if(assmblyName=='hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    BS_genome=BSgenome.Hsapiens.UCSC.hg38
    gen=unique(genome(BS_genome))
  }else if(assmblyName=='t2t'){
    library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
    BS_genome=BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
    seqlevelsStyle(BS_genome) <- "UCSC"
    gen=unique(genome(BS_genome))
  }
  chr=location$chr

  tracksList=list()
  itrack <- IdeogramTrack(genome = gen, chromosome = chr )
  tracksList[[1]]=itrack
  gtrack <- GenomeAxisTrack()
  tracksList[[2]]=gtrack

  grTracks.List=purrr::map(grlist, AnnotationTrack )
  for(x in seq_along(names(grlist))) {grTracks.List[[x]]@name = names(grlist)[x]}

  tracksList=append(tracksList,grTracks.List)

  data(geneModels)
  grtrack <- GeneRegionTrack(geneModels, genome = gen,chromosome = chr, name = "Gene Model", transcriptAnnotation = "symbol", background.title = "brown")
  displayPars(grtrack) <- list(background.panel = "#FFFEDB", col = NULL)
  tracksList=append(tracksList,grtrack)

  # DNA seq ##
  strack <- SequenceTrack(Hsapiens, chromosome = chr)
  tracksList=append(tracksList,strack)

  return(tracksList)
}
