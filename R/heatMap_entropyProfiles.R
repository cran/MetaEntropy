#' Graphically compare entropy profiles
#'
#' This function prints heatmaps of multiple entropy profiles
#' (\code{entropyProfile} objects).
#'
#'
#' @param ... Unquoted \code{entropyProfile} data objects to be compared.
#' @param keepInvariant Logical; whether to retain positions where entropy is
#'   zero. Defaults to \code{FALSE}.
#' @param showAllPeptides Logical; whether the heatmap should show all peptides
#'   encoded by the virus, regardless of whether they have undergone mutations.
#'   Defaults to \code{TRUE}.
#' 
#' @return A \code{ggplot} object representing the heatmap.
#'
#' @importFrom dplyr summarize group_by
#' @importFrom grDevices rgb
#' @importFrom rlang enquos
#'
#' @examples
#' firstWave <- wWater[ wWater$wave == "first", ]
#' thirdWave <- wWater[ wWater$wave == "third", ]
#' # Filter out minor variants whose genotype matches that of the reference genome
#' firstWave <- firstWave[ firstWave$alt_aa_freq <= 0.97, ]
#' thirdWave <- thirdWave[ thirdWave$alt_aa_freq <= 0.97, ]
#' ancestral <- getEntropySignature(firstWave)
#' omicron <-  getEntropySignature(thirdWave)
#' 
#' # High entropies sustained over time in the spike protein, and signals
#' # compatible with pervasive negative selection in Omicron sublineages:
#' heatmap_entropyProfiles(ancestral, omicron)
#' 
#' @export
#
heatmap_entropyProfiles <- function(..., keepInvariant = TRUE, showAllPeptides = TRUE){
	#
	# Capture arguments passed into the ellipsis
	quos <- rlang::enquos(...)
	quos_names <- names(quos) # Capture the argument/list names
	# Validate
	valid_profiles <- list()
	for (i in seq_along(quos)) {
		# Determine profile names in a way compatible with lists,
		# do.call(), and !!!: use the list/argument name if it exists,
		# otherwise fall back to the expression label (for manual inputs)
		if (!is.null(quos_names) && quos_names[i] != "") {
			# For heatmap_entropyProfiles(!!!profiles),
			# quos_names[i] will cleanly pull, e.g. "c1_heart",
			# "c5_liver", etc. (see "intraHostVariants" vignette),
			# directly from the "profiles" list keys.
			expr_name <- quos_names[i]
		} else {
			# For manual inputs like heatmap_entropyProfiles(c1_heart),
			# quos_names will be empty, so the function uses
			# as_label() tu pull "c1_heart".
			expr_name <- rlang::as_label(quos[[i]])
		}
		# Evaluate the expression
		current_obj <- tryCatch(rlang::eval_tidy(quos[[i]]), error = function(e) NULL)
		if (is.null(current_obj)) {
			warning(paste0("Object '", expr_name, "' not found. Skipping."))
			next
		}
		# Perform class check
		obj_class <- class(current_obj)[1]
		if (obj_class == "entropyProfile") {
			valid_profiles[[expr_name]] <- current_obj
		} else {
			warning(paste0("Object '", expr_name, "' is a ", obj_class, ", not entropyProfile. Skipping."))
		}
	}
	#
	# Extract and format data
	profiles <- data.frame()
	for(profile in names(valid_profiles)){
		thisProfile <- valid_profiles[[ profile ]]$Entropy |>
			dplyr::group_by(protein) |>
			dplyr::summarize(entropy = mean(entropy), .groups = "drop") |>
			as.data.frame()
		thisProfile[["stratum"]] <- profile
		profiles <- rbind(profiles, thisProfile)
	}
	# Process
	# Check if zero entropy positions are to be displayed
	if (!keepInvariant){
		profiles <- profiles[profiles$entropy > 0, ]
		minColor <- grDevices::rgb(1, 0.95, 0.95) # very faint red (lowest entropy)
		maxColor = grDevices::rgb(1, 0, 0)
	}
	else{
		minColor <- grDevices::rgb(1, 1, 1) # white (entropy = 0)
		maxColor = grDevices::rgb(1, 0, 0)
	}
	# See if all peptides are to be represented in the heatmap and ensure
	# they are plotted in genome order:
	# Extract ALL proteins from ALL valid profiles. Note that this
	# allows comparison of profiles generated using different ways of
	# encoding the genome (e.g., one with nsp12a and nsp12b, and another
	# with nsp12; see "intraHostVariants" vignette)
	all_reference_peptides <- unique(unlist(lapply(valid_profiles, function(x) x$Genome$CDS$protein)))
	if(showAllPeptides){
		# Use the combined reference list safely
		maturePeptides <- all_reference_peptides
		profiles$protein <- factor(profiles$protein,
					   levels = maturePeptides
		)
	}
	else{
		# Intersect with the combined reference list safely
		maturePeptides <- intersect(all_reference_peptides, profiles$protein)
		profiles$protein <- factor(profiles$protein,
					   levels = maturePeptides)
	}
	#
	# Plot
	return(
	      ggplot2::ggplot(profiles, ggplot2::aes(x = protein, y = stratum, fill = entropy)) +
		      ggplot2::geom_tile(color = "black", linewidth = 0.5) +
		      ggplot2::scale_fill_continuous(low = minColor, high = maxColor, name = "Mean\nEntropy") +
		      ggplot2::scale_x_discrete(drop = FALSE) +
		      ggplot2::theme_minimal() +
		      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
	)
	#
}
