#' Print method for  \code{profileSummary} objects
#'
#' This function formats and prints compact entropy profile summaries
#' (\code{profileSummary} objects), on the console.
#'
#' @param x An object of class \code{profileSummary} created by
#'          \code{\link{summary.entropyProfile}}.
#' @param ... Additional arguments passed to the function.
#'
#' @method print profileSummary
#' 
#' @return Invisibly returns \code{NULL}. This function is used for its side
#'         effect.
#' 
#' @export
#
print.profileSummary <- function(x, ...) {
	gsub("Table: ",
	     "",
	     knitr::kable(x$summary,
			  col.names = c("Feature",
					"Min.",
					"1st Qu.",
					"Median",
					"Mean",
					"3rd Qu.",
					"Max.",
					"SNVs/Kb",
					"non-Syn/Kb"
					),
			  caption = paste0("\n'",
					   paste0(toupper(substr(x$metagenome, 1, 1)),
						  tolower(substr(x$metagenome, 2, nchar(x$metagenome)))
						  ),
					   "' entropy signature.\n\nSummary of entropy and mutations by protein:"),
			  align = "c",
			  digits = 2
	     )

	) |> cat(sep = "\n")
	#
	invisible(NULL)
}
