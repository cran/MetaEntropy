#' Summarize entropy profile
#'
#' Prints a report about an entropy profile (an object of class
#' "entropyProfile").
#'
#' @param object An object of class \code{entropyProfile}.
#' @param ...    Other parameters passed to the function.
#'
#' @rdname summary.entropyProfile
#'
#' @method summary entropyProfile
#' 
#' @return An object of class \code{c("profileSummary", "list")} summarizing
#'         an entropy profile. Intended to be displayed via
#'         \code{print.profileSummary}.
#' 
#' @export
#
summary.entropyProfile <- function(object, ...){
	#
	feature <- c("allCDS", as.character(object$Mutations$protein))
	df <- data.frame(feature	 = feature,
			 min		= numeric(length(feature)),
			 firstQu	= numeric(length(feature)),
			 median		= numeric(length(feature)),
			 mean		= numeric(length(feature)),
			 thirdQu	= numeric(length(feature)),
			 max		= numeric(length(feature)),
			 SNV_Kb		= numeric(length(feature)),
			 nonSyn_Kb	= numeric(length(feature))
	)
	#
	# generic descriptive statistics
	for(df_row in 1:dim(df)[1]){
		if(df$feature[df_row] == "allCDS"){
			df[df_row, 2:7] <- as.numeric(summary(object$Entropy$entropy))
		}
		else{
			df[df_row, 2:7] <- as.numeric(summary(object$Entropy$entropy[object$Entropy$protein == df$feature[df_row]]))
		}
	}
	#
	# Sumarize mutations
	# total (all CDS; syn[1]) and per-protein synonymous mutations
	syn <- c(sum(object$Mutations$syn), object$Mutations$syn)
	# the same for non-synonymous mutations
	nonSyn <- c(sum(object$Mutations$nonSyn), object$Mutations$nonSyn)
	SNVs <- syn + nonSyn
	# All CDS and per-protein lengths
	cdsLengths <- c(sum(object$Mutations$cdsLength), object$Mutations$cdsLength)
	# mutational densities (all CDSs plus every CDS)
	norm_SNVs <- 1000 * SNVs/cdsLengths
	# non-synonymous mutational densities
	norm_nonSyn <- 1000 * nonSyn/cdsLengths
	df$SNV_Kb <- norm_SNVs
	df$nonSyn_Kb <- norm_nonSyn
	#
	lst <- list(summary = df,
			metagenome = deparse(substitute(object)),
			SNVnum = sum(SNVs),
			synNum = sum(syn),
			nonSynNum = sum(nonSyn)
	)
	#
	class(lst) <- c("profileSummary", class(lst))
	#
	return(lst)
	#
}
