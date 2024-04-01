#' @title  Liftover genomic range  (Human genome)
#'
#' @description  Liftover genomic range Human genome and will also add  assembly information to the output genomic ranges. Available for hg19, hg38, t2t genome assemblies. (filter out non-standard chromosomes)
#'
#' @param liftFromTo
#' @param gr.f A genomic range.
#'
#' @return Liftovered GenomicRange object.
#' @export
#' @examples
#' data(hg_38_gr)
#' outGr_h19=liftoverGR(gr.f=hg_38_gr, liftFromTo='hg38.To.hg19')
#' outGr_h19
#'
#'
liftoverGR=function(gr.f, liftFromTo){
  seqlevelsStyle(gr.f) = "UCSC"  # necessary

  if(liftFromTo=='hg38.To.hg19'){
    path = system.file( "extdata", "hg38ToHg19.over.chain", package = "CustomBioInfoFuctionsHumanGenome", mustWork = TRUE) #library(rtracklayer); library(liftOver)
  }else if(liftFromTo=='hg38.To.t2t'){
    path = system.file("extdata", "hg38-chm13v2.over.chain", package = "CustomBioInfoFuctionsHumanGenome", mustWork = TRUE) #library(rtracklayer); library(liftOver)
  }else if(liftFromTo=='hg19.To.hg38'){
    path = system.file("extdata", "hg19ToHg38.over.chain",package = "CustomBioInfoFuctionsHumanGenome", mustWork = TRUE) #library(rtracklayer); library(liftOver)
  }else if(liftFromTo=='hg19.To.t2t'){
    path = system.file("extdata", "hg19-chm13v2.over.chain", package = "CustomBioInfoFuctionsHumanGenome", mustWork = TRUE) #library(rtracklayer); library(liftOver)
  }else if(liftFromTo=='t2t.To.hg19'){
    path = system.file("extdata", "chm13v2-hg19.over.chain", package = "CustomBioInfoFuctionsHumanGenome", mustWork = TRUE) #library(rtracklayer); library(liftOver)
  }else if(liftFromTo=='t2t.To.hg38'){
    path = system.file( "extdata", "chm13v2-hg38.over.chain",package = "CustomBioInfoFuctionsHumanGenome", mustWork = TRUE) #library(rtracklayer); library(liftOver)
  }
  ch = import.chain(path)
  gr.f = unlist(liftOver(gr.f, ch))
  gr.f=addGrInformation(gr.f,strsplit(liftFromTo, split="[.]")[[1]][3])
  return(gr.f)
}
