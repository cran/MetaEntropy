#' Calculates aminoacidic entropy resulting from SNVs present in a specific locus
#'
#' The function is used internally by \code{\link{getEntropySignature}}.
#' It calculates Shannon entropy from categories and frequencies passed to the
#' function.
#'
#' @param variantPosition A list, created  and passed to the function by
#'                        \code{fillPosition}. It contains information on the
#'                        frequencies of each amino acid category observed in a
#'                        virome, as a result of mutations at a given genomic
#'                        position.
#'
#' @return A \code{numeric} value, corresponding to the entropy associated to
#'         the SNVs present at a specific locus.
#'
#' @seealso \code{\link{getEntropySignature}}.
#'
getPosEntropy <- function(variantPosition){
	frecuencias <- variantPosition$frequencies[variantPosition$frequencies > 0]
	entropia <- 0
	for(thisFrequency in frecuencias){
		       		entropia <- entropia + thisFrequency*log(thisFrequency,
										base = length(variantPosition$frequencies))# scale to [0, 1] range.
	       		   }
	return(-entropia)
}
