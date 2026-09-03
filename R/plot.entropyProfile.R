#' Plot entropy signatures
#'
#' Creates entropy charts along a genome.
#'
#' @param x  Object of class \code{entropyProfile}.
#' @param chartType Whether to graph per-protein summaries ("bp"), per-protein
#'                  stripcharts ("stripchart" / "swarm"), or position-wise
#'                  entropy ("entroScan").
#' @param ... Additional arguments passed to the function.
#'
#' @return Unrendered \code{gg/ggplot} object produced by \code{ggplot2}. This
#'                    function is primarily called for its side effects.
#'
#' @importFrom patchwork wrap_plots
#'
#' @examples
#' ancestral <- getEntropySignature(wWater[wWater$wave == "first", ])
#' omicron <- getEntropySignature(wWater[wWater$wave == "third", ])
#' 
#' # Enhanced Spike entropy plus pervasive negative selection in Omicron
#' # sublineages
#' anc_plot <- plot(ancestral, chartType = "stripchart")
#' omi_plot <- plot(omicron, chartType = "stripchart")
#' patchwork::wrap_plots(anc_plot/omi_plot)
#'
#'
#' @export
#'
#' @method plot entropyProfile
#'
plot.entropyProfile <- function(x, chartType = "bp", ...){
	switch(chartType,
	       bp = {ggplot2::ggplot(as.data.frame(x), ggplot2::aes(x = protein, y = entropy)) +
		            ggplot2::geom_boxplot() +
			    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
	       },
	       stripchart = {ggplot2::ggplot(as.data.frame(x), ggplot2::aes(x = protein, y = entropy, col = protein)) +
		       ggplot2::geom_jitter(width = 0.2) +
		       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

	       },
	       entroScan = {ggplot2::ggplot(as.data.frame(x), ggplot2::aes(x = position, y = entropy, col = protein)) +
		       ggplot2::geom_point() +
		       ggplot2::labs(x = "genomic position") +
		       ggplot2::coord_cartesian(xlim = c(1, x$Genome$Length)) +
		       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
	       },
	       swarm = {ggplot2::ggplot(as.data.frame(x), ggplot2::aes(x = protein, y = entropy, col = protein)) +
		       ggbeeswarm::geom_beeswarm() +
		       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
	       }
	)
	
}
