#' @title  Create GenomicRanges object of given the bin size  for human genome.
#'
#' @description Create GenomicRanges object of given the bin size  for human genome.
#' @param assmblyName hg19 or hg38 or t2t
#' @param binSize  size of the genomic block
#' @param chrName  name of the chromosome ((UCSC format)), if not given then all chromosomes are considered (optional)
#' @return GenomicRanges object of given the bin size
#' @export
#' @examples
#' hg_19_Bins.gr=getGbins(assmblyName='hg19', binSize=1000 )
#' hg_19_Bins.gr
#'
#'
getGbins=function(assmblyName, binSize=1000, chrName=chrName){
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

  if(missing(chrName)){
    fgbins=genomeBlocks(BS_genome, chrs = seqnames(BS_genome), width = binSize) # library("Repitools")
  }else{
    fgbins=genomeBlocks(BS_genome, chrs = chrName, width = binSize) # library("Repitools")
  }
  fgbins=addGrInformation(fgbins,assmblyName)
  fgBinSeqs=getSeq(x=BS_genome, fgbins)
  fgbins$CpG_counts=2*oligonucleotideFrequency(DNAStringSet(fgBinSeqs), width=2, step=1)[,'CG']
  return(fgbins)


}
