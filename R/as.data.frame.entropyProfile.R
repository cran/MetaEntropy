#' Coerce \code{entropyProfile} to a Data Frame
#'
#' Function to extract summary information from an \code{entropyProfile}
#' object. This function is internally used for plotting.
#'
#' @param x An object of class \code{entropyProfile}.
#' @param row.names Please see \code{\link{as.data.frame}}.
#' @param optional Please see \code{\link{as.data.frame}}.
#' @param ... Additional arguments passed to the function.
#'
#' @return A \code{data frame} with tabular information on an entropy profile.
#'         This information includes the name of the proteins presenting
#'         mutations, the corresponding genomic positions, and the resulting
#'         entropies in the metagenome.
#'
#' @method as.data.frame entropyProfile
#'
#' @export
#'
as.data.frame.entropyProfile <- function(x, row.names = NULL, optional = FALSE, ...){
	return(data.frame(protein = x$Entropy["protein"],
			  position = x$Entropy["position"],
			  entropy = x$Entropy["entropy"]
               )
	)
}
