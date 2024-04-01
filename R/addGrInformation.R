#' @title  add assembly information to the genomic range
#'
#' @description  This function will add genomic length and assembly name to given genomic ranges (Human genome only). It will remove the non-standard chromosomes from genomic ranges and report the bad genomic ranges for the selected genome assembly.
#' @param gr.f A genomic range.
#' @param assmblyName human genome assembly name i.e., hg19 or hg38 or t2t.
#' @return Input genomic ranges with added assembly information.
#' @export
#' @examples
#'
#' gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),  ranges = IRanges(start = c(1,110,105), end = c(100, 120, 150 )))
#' outGr=addGrInformation(gr.f=gr, assmblyName='hg19')
#'
#'
#'
addGrInformation=function(gr.f, assmblyName){

  if(seqlevelsStyle(gr.f)[1]=='NCBI'){
    seqlevelsStyle(gr.f) <- "UCSC"
    print('Changed seqLevels')
  }


  if(assmblyName=='hg19'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    BS_genome=BSgenome.Hsapiens.UCSC.hg19
    genome(gr.f) = unique(genome(BSgenome.Hsapiens.UCSC.hg19))
  }else if(assmblyName=='hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    BS_genome=BSgenome.Hsapiens.UCSC.hg38
    genome(gr.f) = unique(genome(BSgenome.Hsapiens.UCSC.hg38))
  }else if(assmblyName=='t2t'){
    library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
    BS_genome=BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
    seqlevelsStyle(BS_genome) <- "UCSC"
    #genome(gr.f) = unique(genome(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0))

  }

  gr.f=keepStandardChromosomes(gr.f, pruning.mode="coarse")
  seqlevels(gr.f, pruning.mode="coarse")=seqlevels(BS_genome)[seq(1,24)]
  seqlengths(gr.f)=seqlengths(BS_genome)[seq(1,24)]
  isCircular(gr.f)=isCircular(BS_genome)[seq(1,24)]

  BadGr=which(end(gr.f) > seqlengths(BS_genome)[as.character(seqnames(gr.f))])
  if(length(BadGr)>0){
    print('bad GRs')
    print(gr.f[BadGr])
    stop("This is an error message: please check assembly of GR")
  }

  return(gr.f)
}
