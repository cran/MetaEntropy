#' Create (empty) object of class "entropyProfile"
#'
#' This function is intended primarily for internal use by
#' \code{\link{getEntropySignature}}.
#'
#' @param polymorphisms  A data frame. Please see Details and Examples in
#'        documentation for \code{\link{getEntropySignature}}.
#' @param position Name of the \code{polymorphisms}'s column that indicates SNV
#'                 locations in the genome.
#' @param linkage Information on linked positions.
#' @param ref Column name with reference bases.
#' @param alt Column name with the alternative bases observed in the
#'            metagenome.
#' @param protein Name of the column carrying protein names.
#' @param aa_position Name of the column that indicates the protein positions
#'                    of the mutated amino acids.
#' @param ref_aa Name of the column that carries the reference amino acids.
#' @param alt_aa Name of the column carrying alternative amino acids observed
#'               in the metagenome.
#' @param alt_aa_freq Name of the column giving the frequencies of alternative
#'                 amino acids in the metagenome.
#' @param entropies \code{NA_REAL_} (double numeric/real vector to hold entropy
#'                  values).
#' @param genome A list providing CDS data and length of the reference genome.
#'
#' @details
#' The documentation for \code{\link{getEntropySignature}} details the type of
#' input needed to create a profile. \code{entropyProfile} uses the same parameters as
#' \code{getEntropySignature}, with the exception of \code{categories} and
#' \code{entropies}.
#'
#' @return An (empty) object of class \code{entropyProfile}.
#'
#' @seealso \code{\link{getEntropySignature}}.
#'
#' @export
#
entropyProfile <- function(polymorphisms,
			   position = "position",
			   linkage = "linkage",
			   ref = "ref",
			   alt = "alt",
			   protein = "protein",
			   aa_position = "aa_position",
			   ref_aa = "ref_aa",
			   alt_aa = "alt_aa",
			   alt_aa_freq = "alt_aa_freq",
			   entropies = NA_real_,
			   genome = mn908947.3){
	#
	# preserve original data (some point mutions will not be mirrored in
	# the entropy profile due to linkage)
	snvs <- polymorphisms
	# standardize labels:
	names(snvs)[ names(snvs) == position ] <- "position"
	names(snvs)[ names(snvs) == linkage ] <- "linkage"
	names(snvs)[ names(snvs) == ref ] <- "ref"
	names(snvs)[ names(snvs) == alt ] <- "alt"
	names(snvs)[ names(snvs) == protein ] <- "protein"
	names(snvs)[ names(snvs) == aa_position ] <- "aa_position"
	names(snvs)[ names(snvs) == ref_aa ] <- "ref_aa"
	names(snvs)[ names(snvs) == alt_aa ] <- "alt_aa"
	names(snvs)[ names(snvs) == alt_aa_freq ] <- "alt_aa_freq"
	#
	# extract coding DNA sequence data (simplify code below)
	cds <- genome$CDS
	#
	# account for linked positions
	linked <- FALSE
	for(snv in 1:dim(polymorphisms)[1]){
		if(is.na(polymorphisms[,linkage][snv])){
		# nothing to do
		}
		else{
			# what positions are linked?
			haplotype_start <- polymorphisms[,position][snv]
			inside <- TRUE
			while(inside){
				haplotype_end <- polymorphisms[ ,linkage][snv]	# e.g.: snv 10 in position 234 linked to snv 11 in position 235
										#    	wave,	position,	linkage
										#     9	xxx,	100,		NA
										#    10	xxx,	234,		235
										#    11	xxx,	235,		234
										#    12	xxx,	583,		NA
										# haplotype_end = x[,linkage][10] = 235
				snv <- snv + 1# move one snv downstream (e.g., snv 11)
				if(polymorphisms[ ,linkage][snv] < polymorphisms[ ,position ][snv]){# Linked to upstream position only
					inside <- FALSE
					snv <- snv + 1# move one snv downstream (e.g. snv 12)
				}
				else{
					# Do nothing (still inside the haplotype)
				}
			}
			# Colapse linked positions; map residue substitution to haplotype_start
			if(haplotype_end - haplotype_start == 1){# 2 codon positions linked
				polymorphisms <- polymorphisms[ -(snv-1),]# e.g. drop row 11
			}
			else{# All 3 positions linked
				polymorphisms <- polymorphisms[ -c(snv-1, snv-2),]
			}
		}
	}
	#
	# proteins with records
	proteins <- character()
	for(posicion in unique(polymorphisms[ ,position])){
		proteins <- c(proteins, unique(polymorphisms[ ,protein][polymorphisms[ ,position] == posicion]))
		# unique() is to account for mutiple haplotypes.
		# E.g.:
		# > polymorphisms[ polymorphisms[,position] == 23429 ,]
		#     wave position linkage ref alt protein aa_position ref_aa alt_aa alt_aa_freq
		# 50 first    23429      NA   G   T       S         623      A      S     0.09884
		# 51 first    23429      NA   G   A       S         623      A      T     0.12706
		# These two rows make a joint contribution to the entropy at
		# position 23429, so the profile will contain a single row for
		# position 23429 (ergo a unique protein label, "S").
	}
	# keep order of proteins (most likely genomic order)
	proteins <- factor(proteins, levels = unique(proteins[!is.na(proteins)]))
	#
	#
	# Structure that carries the profile
	perfil <- list(SNVs = snvs,
		       Entropy = data.frame(position = unique(polymorphisms[,position, drop = T]),# again account for multi-haplotype polymorphisms
					    protein = proteins,
					    entropy = entropies
		       ),
		       Mutations = data.frame(protein = unique(proteins), # proteins with records
					      cdsLength = cds[ which(cds[ , "protein"] %in% unique(proteins)), "end"] - cds[ which(cds[ , "protein"] %in% unique(proteins)), "start"] - 1,
					      syn = numeric(length = length(unique(proteins))),
					      nonSyn = numeric(length = length(unique(proteins)))
					      ),
		       Genome = genome 
	)
	#
	# Summarize mutation data. ("polymorphisms" has been already processed above for linked haplotypes)
	for(cdsRegion in unique(perfil$Entropy$protein)){
		perfil$Mutations$syn[perfil$Mutations$protein == cdsRegion] <- sum(polymorphisms[ polymorphisms$protein == cdsRegion, "ref_aa"] == polymorphisms[ polymorphisms$protein == cdsRegion, "alt_aa"])
		perfil$Mutations$nonSyn[perfil$Mutations$protein == cdsRegion] <- sum(polymorphisms[ polymorphisms$protein == cdsRegion, "ref_aa"] != polymorphisms[ polymorphisms$protein == cdsRegion, "alt_aa"])
	}
	#
	#
	class(perfil) <- c("entropyProfile", class(perfil))
	return(perfil)
}
