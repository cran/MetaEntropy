#' Summarize mutations
#'
#' Displays SNVs, and corresponding protein mutations, at specific genomic
#' positions.
#'
#' @param profile  An object of class \code{entropyProfile}.
#' @param positions A vector with genome positions relative to the reference
#'        genome.
#'
#' @return An object of class \code{c("tidyMutations", "data.frame")},
#'         containing summary information about user-supplied genomic
#'         positions. This information includes the mutations themselves
#'         relative to the reference genome, their positions within it, and the
#'         corresponding abundances in the virome. Intended to be displayed by
#'         \code{print.tidyMutations}.
#'
#' @details The user provides a list of genome positions and the function
#'          prints the mutations associated with them.
#'          The output format is "ref_res###alt_res /
#'          protein:ref_res###alt_res", where ref_res is the residue (eiter
#'          nucleotide or aminoacid) in the reference strain, alt_res is the
#'          alternative residue in the metagenome, "###" is the position
#'          (either nucleotide or aminoacid) where the mutation was observed,
#'          and "protein" is the name of the affected protein.
#'
#' @seealso \code{\link{getEntropySignature}}.
#'
#' @examples
#' 
#' # High entropy at the RBD in Omicron lineages
#' omicron <- getEntropySignature(wWater[wWater$wave == "third", ])
#' plot(omicron, chartType="stripchart")
#'
#' # Identify the high-entropy positions
#' omicron$Entropy$position[ omicron$Entropy$entropy > 0.3 ]
#' #[1] 22882 22898 22917 23013 23040 23048 23055 23063
#' 
#' # Get a descriptive table
#' showMutations(omicron, c(22882, 22898, 22917, 23013, 23040, 23048, 23055, 23063))
#'
#'
#' @export
#
showMutations <- function(profile, positions){
	#
	if(!(inherits(profile, "entropyProfile"))){
		stop("The \"profile\" parameter must be an object of class \"entropyProfile\"")
	}
	if(!is.numeric(positions) | !is.vector(positions)){
		stop("The \"positions\" parameter must be a numeric vector (e.g. c(22882, 22898, 22917))")
	}
	#
	# See if "positions" have mutations
	if(any(positions %in% profile$SNVs$position)){ # at least one position displays mutations:
		# See if there are invariant positions also
		if(length(positions)# number of positions passed to the function
		   >
		   sum(positions %in% profile$SNVs$position)){# number of such positions having mutations
			withoutRecords <- positions[ !(positions %in% profile$SNVs$position) ]
			message(paste("No mutations recorded at:", paste(withoutRecords, collapse = ", "), "\n"))
			positions <- positions[ positions %in% profile$SNVs$position ]
		}
	}
	else{
		message(paste("No mutations recorded at:", paste(positions, collapse = ", "), "; returning NULL"))
		return(NULL)
	}
	#
	# Account for multi-haplotype positions
	snv_rows <- which(profile$SNVs$position %in% positions)
	#  
	#
	df <- data.frame(snv = paste0(profile$SNVs$ref[ snv_rows ],
					   profile$SNVs$position[ snv_rows ],
					   profile$SNVs$alt[ snv_rows ]
					   ),
			 phenotype = paste0(profile$SNVs$protein[ snv_rows ],
					    ":",
					    profile$SNVs$ref_aa[ snv_rows ],
					    profile$SNVs$aa_position[ snv_rows ],
					    profile$SNVs$alt_aa[ snv_rows ]
					    ),
			 abundance = round(profile$SNVs$alt_aa_freq[ snv_rows ] * 100)
	)
	#
	# Set class for pretty console output ("print.tydyMutations.R")
	class(df) <- c("tidyMutations", class(df))
	#
	return(df)
	#
}
