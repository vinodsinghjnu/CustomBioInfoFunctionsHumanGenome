#' @title  Get genomic range of ProteinCodingGenes promoters (upst to -200)
#'
#' @description  This will return a promoters / genes genomic ranges for a given human assembly.
#' @param upst Upstream of the promoter. Promoter region will be "upst" to -200 to TSS. DEFAULT value is 500
#' @param assmblyName human genome assembly name i.e., hg19 or hg38
#' @return  genomic ranges for protein coding genes and their promoter for a give assembly.
#' @export
#' @examples
#'
#'
#' outGr=getPromGRange(upst=500, assmblyName='hg19')
#' outGr$proteinCoding.proms
#' outGr$proteinCoding.genes
#'
#'
#'

getPromGRange=function(upst=500, assmblyName){
   if(assmblyName=='hg19'){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene #shorthand (for convenience)
   }else if(assmblyName=='hg38'){
     library(TxDb.Hsapiens.UCSC.hg38.knownGene)
     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene #shorthand (for convenience)
   }
    prom_gr=addGrInformation(promoters(txdb, upstream=upst, downstream=0), assmblyName)
    #genes_gr=addGrInformation(unlist(sort(genes(txdb))), assmblyName)
    genes_gr=addGrInformation(sort(unlist(genes(txdb, single.strand.genes.only=FALSE))), assmblyName)
    genes_gr$gene_id=names(genes_gr)



    #library("org.Hs.eg.db")
    GeneMapID.df=select(org.Hs.eg.db, keys=genes_gr$gene_id, keytype = 'ENTREZID',
                        columns = c("ENTREZID","SYMBOL","GENETYPE","GENENAME", 'ENSEMBL'))
    GeneMapID.df=GeneMapID.df[! duplicated(GeneMapID.df$ENTREZID),]
    rownames(GeneMapID.df)=GeneMapID.df$ENTREZID
    genes_gr$GenesMap= GeneMapID.df[genes_gr$gene_id,]

    genes_Promoter_gr=promoters(genes_gr, upstream = upst)
    ProteinCoding_genes_Promoter_gr=genes_Promoter_gr[which(genes_Promoter_gr$GenesMap$GENETYPE=="protein-coding"),]
    ProteinCodingGenes_UCSC_knownGenes_detailed=genes_gr[which(genes_gr$GenesMap$GENETYPE=="protein-coding"),]
    return(list(proteinCoding.proms=ProteinCoding_genes_Promoter_gr, proteinCoding.genes=ProteinCodingGenes_UCSC_knownGenes_detailed))
}




