#' Summarize variants and frequencies at a genome position
#'
#' This function is used internally by \code{getEntropySignature()}.
#' It creates a vector (\code{aminoAcids}) listing the amino acids
#' observed in a virome at a particular position under analysis, including the
#' reference amino acid, another vector (\code{frequencies}) with the
#' corresponding frequencies, and returns them combined in a data frame.
#'
#' @param variants A data frame, similar to the \code{polymorphisms} argument
#'                 of \code{\link{getEntropySignature}}, but containing
#'                 information on a single genome position.
#' @param ref_aa Name of the column that carries reference amino acids.
#' @param alt_aa Name of the column carrying alternative amino acids observed
#'               in the metagenome.
#' @param alt_aa_freq Name of the column giving the frequencies of alternative
#'                 amino acids.
#'
#' @return A \code{data frame} describing the variability (different amino
#'         acids its frequencies) observed at a specific locus.
#'
#' @seealso \code{\link{getEntropySignature}}.
#'
#' @keywords internal
#
createPositionSummary <- function(variants, ref_aa, alt_aa, alt_aa_freq){
	aminoAcids <- character(length = dim(variants)[1] + 1)
	aminoAcids[1] <- variants[,ref_aa][1]
	frequencies <- numeric(length = dim(variants)[1] + 1)
	for(variante in 1:length(variants[,alt_aa])){
		aminoAcids[variante + 1] <- variants[,alt_aa][variante]
		frequencies[variante + 1] <- variants[,alt_aa_freq][variante]
	}
	frequencies[1] <- 1 - sum(frequencies[frequencies != 0])
	return(data.frame(aminoAcids = aminoAcids, frequencies = frequencies))
}
