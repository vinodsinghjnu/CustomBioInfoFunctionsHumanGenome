#' @title  Generate all possible DNA sequences of a \href{https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide}{Ambiguous nucleotide sequence}
#'
#' @description Generate all possible DNA sequences of a \href{https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide}{Ambiguous nucleotide sequence}
#' @param pat Ambiguous nucleotide sequence
#' @return A vector of all possible DNA sequences for a given Ambiguous nucleotide sequence
#' @export
#' @examples
#' DNA_seqs=DNASeqsForPattern(pat='NYYN')
#' DNA_seqs
#'
#' @seealso code{\link[DECIPHER]{Disambiguate}}
#'
#'
DNASeqsForPattern=function(pat){
  possibleCombination=as.vector(unlist(Disambiguate(DNAStringSet(pat))))
  return(possibleCombination)
}
