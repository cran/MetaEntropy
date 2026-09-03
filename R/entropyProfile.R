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
	# Account for ORF superposition (positions are multiplexed internally)
	# Implements "multiplexed" positions by adding decimals, e.g. position
	# 100, becomes 100, 100.1, 100.2, etc. ("pseudo mapping").
	#
	# A "general" example, with "linkage" (not a real example);
	# All three positions in a codon display mutations
	#
	#       position linkage ref alt  protein   codon    ref_aa alt_aa  freq
	# 50    1        2       G   T    X         123      A      B       0.09884
	# 51    1        2       A   A    Y         456      M      N       0.09884
	# 52    2        3       G   T    X         123      A      B       0.09884
	# 53    2        3       C   G    Y         456      M      N       0.12706
	# 54    3        2       A   A    X         123      A      B       0.09884
	# 55    3        2       C   G    Y         456      M      N       0.12706
	#
	# Ergo, have to "displace" either all the positions associated with protein X or all
	# the positions associated with protein Y.
	#
	# 1. See which positions are repeated
	are_duplicate <- duplicated(polymorphisms[ ,position])
	# returns a logic pointing to cells 50 to 55.
	# 
	duplicates <- unique(polymorphisms[ ,position][are_duplicate])
	# returns "c(1, 2, 3)"
	#
	# 2. See which positions are  repeated due to overlaps ("protein" column will be labelled differently).
	for(this_duplicate in duplicates){
		# check if a single or multiple proteins are listed for the repeated position 
		affected_proteins <- polymorphisms[ polymorphisms[ ,position] == this_duplicate, protein]
		# when this_duplicate is 1 affected_proteins is "c(X, Y)"
		n_affected_proteins <- length(unique(affected_proteins))
		# when this_duplicate is 1, n_affected_proteins is "2" 
		if( n_affected_proteins > 1 ){# there is an overlap; multiplex it (true when this_position equals 1, 2, and 3)
			# See how many times the position is repeated
			repetitions <- length(polymorphisms[ polymorphisms[,position] == this_duplicate ,position])
			# equals "2" for positions 1, 2 and 3
			# Calculate displacements
			displacement <- 1/repetitions # e.g. 4 repetitions gives 0.25.; then
			                              # "x" becomes "x", "x.25", "x.50", "x.75" (never reaches x + 1)
			# diplacement is 0.5 for positions 1, 2 and 3
			# multiplex
			toMultiplex <- which( polymorphisms[,position] == this_duplicate )# a numerical vector, indexing
			                                                                  # the positions to be multiplexed
			# equals "c(50, 51)" when this_duplicate is "1"
			for(displaceThis in toMultiplex[-1]){# [-1] remove the first position of the vector
				# loop starts at "52" when this position is "1"
				polymorphisms[,position][displaceThis] <- polymorphisms[,position][displaceThis - 1] + displacement
				# Using "[displaceThis - 1]" keeps adding "displacement" the the previous member of the
				# When this_position equals "1", our example becomes:
				#
				#       position   linkage ref alt  protein   codon    ref_aa alt_aa  freq
				# 50    1          2       G   T    X         123      A      B       0.09884
				# 53    1.5        2       A   A    Y         456      M      N       0.09884
				# 51    2          3       G   T    X         123      A      B       0.09884
				# 54    2          3       C   G    Y         456      M      N       0.12706
				# 52    3          2       A   A    X         123      A      B       0.09884
				# 55    3          2       C   G    Y         456      M      N       0.12706
				#
				# So, we have to multiplex the linked position informing on protein "Y" (position[52]),
				# because in the next round, when this_positions will be "2", it will be set to 2.5.
				#
				# Check if linked
				this_linkage <- polymorphisms[ , linkage][displaceThis]
				if(!is.na(this_linkage)){
					polymorphisms[ , linkage][displaceThis] <- this_linkage + displacement
				}
				# When thisposition is 1 displaceThis is "52", so the above results in
				#
				#       position   linkage   ref alt  protein   codon    ref_aa alt_aa  freq
				# 50    1          2         G   T    X         123      A      B       0.09884
				# 51    1.5        2.5       A   A    Y         456      M      N       0.09884
				# 52    2          3         G   T    X         123      A      B       0.09884
				# 53    2          3         C   G    Y         456      M      N       0.12706
				# 54    3          2         A   A    X         123      A      B       0.09884
				# 55    3          2         C   G    Y         456      M      N       0.12706
				#
				# And after processing all 3 positions:
				#
				#       position   linkage   ref alt  protein   codon    ref_aa alt_aa  freq
				# 50    1          2         G   T    X         123      A      B       0.09884
				# 51    1.5        2.5       A   A    Y         456      M      N       0.09884
				# 52    2          3         G   T    X         123      A      B       0.09884
				# 53    2.5        3.5       C   G    Y         456      M      N       0.12706
				# 54    3          2         A   A    X         123      A      B       0.09884
				# 55    3.5        2.5       C   G    Y         456      M      N       0.12706
				#
			}
		}
		else{# no overlap, ergo do nothing
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
		       Entropy = data.frame(position = unique(polymorphisms[,position]),# again account for multi-haplotype polymorphisms
					    protein = proteins,
					    entropy = entropies
		       ),
		       Mutations = data.frame(protein = unique(proteins), # proteins with records
					      cdsLength = cds[ which(cds[ , "protein"] %in% unique(proteins)), "end"] - (cds[ which(cds[ , "protein"] %in% unique(proteins)), "start"] - 1),
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
	#
	return(list(Perfil = perfil,
		    Polymorphisms_multiplexed = polymorphisms
		    )
	)
}
