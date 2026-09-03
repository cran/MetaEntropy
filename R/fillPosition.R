#' Translate amino acid frequencies into category frequencies
#'
#' The function is used internally by \code{\link{getEntropySignature}}.
#' It creates a storage list by \code{createStorage}, and loads on it the
#' frequency of each amino acid category based on the data contained in a data
#' frame passed to the function (\code{positionsSummary} parameter).
#'
#' @param positionSummary A data frame created by \code{\link{createPositionSummary}}.
#' @param categories A character string indicating which category scheme to
#'                   use. Similar to the \code{categories} parameter of
#'                   \code{\link{getEntropySignature}}.
#'
#' @return A list with information on the frequencies of each amino acid
#'         category observed in a virome in a specific locus. The list contains
#'         a \code{character} vector for each amino acid category, and a named
#'         \code{numeric} vector containing the frequency of each category in
#'         the metagenome. 
#'
#' @seealso \code{\link{getEntropySignature}}.
#'
fillPosition <- function(positionSummary, categories){
	aminoAcids <- positionSummary$aminoAcids
	frequencies <- positionSummary$frequencies
	#
	variantPosition <- createStorage(categories)
	#
	for(aminoAcid in 1:length(aminoAcids)){
		for(categoria in names(variantPosition[names(variantPosition) != "frequencies"])){
			if(!is.na(match(aminoAcids[aminoAcid], variantPosition[[categoria]]))){
				variantPosition$frequencies[categoria] <- variantPosition$frequencies[categoria] + frequencies[aminoAcid]
			}
		}
	}
	return(variantPosition)
}
