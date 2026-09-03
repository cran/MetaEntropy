#' Print method for  \code{tidyMutations} objects
#'
#' This function formats and prints compact mutation summaries
#' (\code{tidyMutations} objects) created by showMutations(), on the console.
#'
#' @param x An object of class \code{tidyMutations} created by
#'          \code{\link{showMutations}}.
#' @param ... Additional arguments passed to the function.
#' 
#' @return Invisibly returns \code{NULL}. Called for side effect.
#'
#' @method print tidyMutations
#' 
#' @export
#' 
#'
print.tidyMutations <- function(x, ...) {
	knitr::kable(x,
		     col.names = c(names(x)[1:2], "abundance (%)"),
		     align = "c") |> cat(sep = "\n")
	invisible(NULL)
}
