#' @title  oligo-nucleotide counts in within a genomic context.
#'
#' @description oligo-nucleotide counts in within a genomic context.
#' @param contextGr GenomicRange object of the genomic context within which oligo-nucleotides has to be counted.
#' @param oligoType dinucs or trinucs or tetranucs
#' @param ignore.strand genomic context strand information should be considered. Default: FALSE
#' @param assmblyName human genome assembly name (hg19 or hg38 or t2t). Default: hg19
#' @return a vector of oligo-nucleotide counts
#' @export
#' @examples
#' data(hg_38_gr)
#' oligonucs.Counts=context_oligonucsCounts(contextGr=hg_38_gr, oligoType='trinucs', ignore.strand=FALSE, assmblyName='hg38')
#'
#'
#'
context_oligonucsCounts=function(contextGr, oligoType, ignore.strand=FALSE,  assmblyName='hg19')
{
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

  if(ignore.strand==FALSE)
  {
    ContextGr_pos=contextGr
    strand(ContextGr_pos)='+'
    ContextGr_neg=contextGr
    strand(ContextGr_neg)='-'
  }else if(ignore.strand==TRUE){
    ContextGr_pos=contextGr[which(strand(contextGr)=='+')]
    ContextGr_neg=contextGr[which(strand(contextGr)=='-')]
  }

  ContextGrSeqs_pos=getSeq(x=BS_genome, ContextGr_pos )
  ContextGrSeqs_neg=getSeq(x =BS_genome, ContextGr_neg )

  if(oligoType=='dinucs'){
    oligonucs_Counts_pos=colSums(dinucleotideFrequency(ContextGrSeqs_pos))
    oligonucs_Counts_neg=colSums(dinucleotideFrequency(ContextGrSeqs_neg))
  }else if(oligoType=='trinucs'){
    oligonucs_Counts_pos=colSums(trinucleotideFrequency(ContextGrSeqs_pos))
    oligonucs_Counts_neg=colSums(trinucleotideFrequency(ContextGrSeqs_neg))
  }else if(oligoType=='tetranucs'){
    oligonucs_Counts_pos=colSums(oligonucleotideFrequency(ContextGrSeqs_pos, width=4, step=1))
    oligonucs_Counts_neg=colSums(oligonucleotideFrequency(ContextGrSeqs_neg, width=4, step=1))
  }

  oligonucs_Counts=oligonucs_Counts_pos + oligonucs_Counts_neg[names(oligonucs_Counts_pos)]

  return(oligonucs_Counts)
}
