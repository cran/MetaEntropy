#' Evaluates Entropy Hotspot
#'
#' Graphical and formal analyses of contiguous amino acids.
#'
#' @param profile An object of class \code{entropyProfile}.
#' @param boundaries Numeric vector with the first and last genomic positions
#'        of the region to be evaluated. To be set interactively if not
#'        provided.
#' @param chartType Chart type; either "boxplot", "stripchart" or "swarm".
#'
#' @return \code{htest} object. This function is called primarily for its side
#'          effects.
#'
#' @details The query stretch (*e.g.* a protein domain with neutralizing
#'          epitopes) is compared against the full set of proteins. Hot spot
#'          boundaries should be indicated relative to the reference genome
#'          used in variant calling.
#'
#' @seealso \code{\link{getEntropySignature}}.
#'
#' @examples
#' omicron <- getEntropySignature(wWater[wWater$wave == "third", ])
#' 
#' # Entrpy hotspot at SARS-CoV-2 receptor binding domain
#' assessHotSpot(omicron, c(22517, 23186), chartType = "swarm")
#'
#' @export
#'
assessHotSpot <- function(profile, boundaries, chartType = "boxplot"){
	#
	correctParameters <- c(perfil = FALSE, boundaries = FALSE)
	while(!correctParameters["perfil"] & !correctParameters["boundaries"]){
		# delimit boundaries:
		if(missing(boundaries)){
			while(!correctParameters["boundaries"]){
				# read boundaries
				boundaries <- readline(prompt = "Enter boundaries separated by a coma: ")
				# inspect input
				incorrectInput <- tryCatch(
							   expr = {
								   boundaries <- as.numeric(unlist(strsplit(boundaries, ",")))# warns if non numeric
								   list(boundaries = boundaries, warned = FALSE)
							   }, warning = function(w){
								   list(boundaries = NULL, warned = TRUE)
							      }
				)
				# check if valid
				if(isTRUE(incorrectInput$warned)){
					message("Boundaries must be numbers separated by commas (e.g. \"22518, 23185\").")
				}
				else{
					correctParameters["boundaries"] <- TRUE
				}
			}
		}
		# check profile
		if(!(inherits(profile, "entropyProfile"))){
			stop("The \"profile\" parameter must be an object of class \"entropyProfile\"")
		}
		else{
			correctParameters["perfil"] <- TRUE
		}
	}
	#
	isQuery <- function(positions){
		positions >= boundaries[1] & positions <= boundaries[2]
	}
	is_hot_spot <- factor(length(profile$Entropy$positions), levels = c("no", "yes"))
	is_hot_spot[isQuery(profile$Entropy$position)] <- "yes"
	is_hot_spot[!isQuery(profile$Entropy$position)] <- "no"
	#
	# graphical analysis 
	switch(chartType,
	       boxplot = {
		       graphics::boxplot(profile$Entropy$entropy ~ is_hot_spot,
					 ylab = "entropy",
					 xlab = "hot spot"
		       )
	       },
	       stripchart = {
		       df <- cbind(profile$Entropy, data.frame(hot_spot = is_hot_spot))
		       graphics::stripchart(entropy ~ is_hot_spot, data = df,
					    vertical = TRUE,
					    method = "jitter",
					    pch = 19,
					    xlab = " hot spot"
		       )
	       },
	       swarm ={
		       beeswarm::beeswarm(profile$Entropy$entropy ~ is_hot_spot,
					  pch = 19,
					  ylab = "entropy",
					  xlab = "hot spot"
		       )
	       }
	)
	#
	# formal analysis
	entropy  <- profile$Entropy$entropy
	hot_spot <- is_hot_spot
	stats::kruskal.test(entropy ~ hot_spot)
}
