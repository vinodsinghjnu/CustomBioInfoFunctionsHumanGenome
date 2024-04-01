#' @title  Creates GenomicRanges object of bins of user specified CpG counts for human genome.
#'
#' @description Creates a GenomicRange of genomic blocks of user specified CpG counts. (human genome)
#' @param assmblyName hg19 or hg38 or t2t
#' @param CGs_perBin  CpG counts in a genomic block/bin. (even number)
#' @param addSeq  if sequence of the bin is required (Default: FALSE)
#' @return GenomicRanges object of bins with user specified CpG counts
#' @export
#' @examples
#' hg_19_CpGBins.gr=Gbin_ByCGcnts(CGs_perBin=1000, assmblyName='hg19' )
#' hg_19_CpGBins.gr
#'
#'
Gbin_ByCGcnts =function(CGs_perBin, assmblyName, addSeq=FALSE){

  CGs_perBin=CGs_perBin/2

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

  CG_Loc = sort(vmatchPattern(DNAString("CG"), BS_genome),ignore.strand=TRUE ) # library(Biostrings)
  #print(head(CG_Loc))
  CG_Loc=CG_Loc[!(strand(CG_Loc) == '-')]
  CG_Loc=addGrInformation(CG_Loc,assmblyName)
  #print(head(CG_Loc))
  CG_Loc_grlist = split(CG_Loc, seqnames(CG_Loc))

  retrunBin=function(gr)
  {
    seqBreaks=seq(1,length(gr), by = CGs_perBin)
    seqBreaks=seqBreaks[-length(seqBreaks)]
    binGr=GRanges(seqnames=rep(seqnames(gr)[1],length(seqBreaks)), ranges = IRanges(start=start(gr)[seqBreaks], end=end(gr)[seqBreaks+CGs_perBin-1]), strand=rep('*',length(seqBreaks)) )
  }

  CG_bin_chrLst=GRangesList(lapply(seq(1,24), function(x) retrunBin(CG_Loc_grlist[[x]])))
  CG_bin=unlist(CG_bin_chrLst)
  if(addSeq==TRUE){
    CG_bin$seq=getSeq(x=BS_genome, CG_bin)
    CG_bin$CpG_counts=2*oligonucleotideFrequency(DNAStringSet(CG_bin$seq), width=2, step=1)[,'CG']
  }
  CG_bin=addGrInformation(CG_bin, assmblyName )
  return(CG_bin)

}
