#' Build a structure representing amino acid categories. 
#' 
#' The function is used internally by \code{\link{fillPosition}}, which in
#' turn is an auxiliary function of \code{\link{getEntropySignature}}.
#' It creates a list with one element for each amino acid category, named
#' according to the categories used (e.g., "aliphatic", "aromatic", etc.).
#' Each element contains a set of amino acids identified by one-letter codes.
#' The list also includes an element containing an empty numeric vector, whose
#' names correspond to the labels of each category.
#' This vector is to be populated with the frequency of each category at a
#' given genomic position, by the `fillPosition` function.
#'
#' @param categories A character string. Similar to the \code{categories}
#'                   parameter of \code{\link{getEntropySignature}}.
#'
#' @return A list with one element (\code{character} vector) for amino acid
#'         category and an element (empty named \code{numeric} vector) to be
#'         loaded with the frequency in the metagenome of each amino acid category.
#'
#' @seealso \code{\link{getEntropySignature}}.
#'
createStorage <- function(categories){ # "sensitive"/"robust"
	switch(categories,
	       robust = catsAndFreqs <- list(
						aliphatic = c("A","V","L","I","M","C"),
						aromatic = c("F","W","Y","H"),
						polar = c("S","T","N","Q"),
						positive = c("K","R"),
						negative = c("D","E"),
						special = c("G","P"),
						STOP = "*",
						frequencies = stats::setNames(# create an empty, named object
								       numeric(length = 7),# the object
								       c(
									 "aliphatic", "aromatic", "polar", "positive", "negative", "special", "STOP"# the names
								        )
						                      )
	                                       ),
	       sensitive = catsAndFreqs <- list(
					     phe = "F",
					     trp = "W",
					     tyr = "Y",
					     met = "M",
					     leu = "L",
					     val = "V",
					     ile = "I",
					     arg = "R",
					     lys = "K",
					     gln = "Q",
					     glu = "E",
					     asp = "D",
					     asn = "N",
					     ser = "S",
					     gly = "G",
					     ala = "A",
					     thr = "T",
					     pro = "P",
					     cys = "C",
					     his = "H",
					     STOP = "*",
					     frequencies = stats::setNames(
								    numeric(length = 21),
								    c(
								      "phe", "trp", "tyr", "met", "leu", "val", "ile", "arg", "lys", "gln", "glu", "asp", "asn", "ser", "gly", "ala", "thr", "pro", "cys", "his", "STOP")# the names
						                   )
			                    )
	)
	return(catsAndFreqs)
}
