#' @title  Create empty chromosome GenomicRange object for a given human genome assembly
#'
#' @description Create empty chromosome GenomicRange object for a given human genome assembly (for standard chromosomes )
#' @param assmblyName hg19 or hg38 or t2t
#' @return empty chromosome GenomicRange
#' @export
#' @examples
#' hg_19_Chr.gr=emptyChrGranges('hg19')
#' hg_19_Chr.gr
#'
#'
emptyChrGranges=function(assmblyName){
  if(assmblyName=='hg19'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    BS_genome=BSgenome.Hsapiens.UCSC.hg19
  }else if(assmblyName=='hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    BS_genome=BSgenome.Hsapiens.UCSC.hg38
  }else if(assmblyName=='t2t'){
    library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
    BS_genome=BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
    seqlevelsStyle(BS_genome) <- "UCSC"
  }
  emptyChrGr=GRanges(names(BS_genome), IRanges(start=1, end=seqlengths(BS_genome)), strand = '*' )

  emptyChrGr=addGrInformation(emptyChrGr,assmblyName)

  chrSeqs=getSeq(BS_genome, emptyChrGr)
  names(chrSeqs)=seqnames(emptyChrGr)
  emptyChrGr$Seqs=chrSeqs
  return(emptyChrGr)
}

