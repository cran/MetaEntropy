## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----loadPackages, results = "hide"-------------------------------------------
lapply(c("MetaEntropy", "ggplot2", "patchwork", "tidyr"),
       library, character.only = TRUE
)

## ----computeSNVs--------------------------------------------------------------
addmargins(
	   table(intraHostVariants$case,
		 intraHostVariants$organ
		 ),
	   margin = 2
)

## ----widenData, results = "hide"----------------------------------------------
wide_ihv <- tidyr::pivot_wider(intraHostVariants,
		 names_from = organ,
		 values_from = alt_aa_freq,
		 # Use 0 (instead of NA) for undetected
		 # mutations.
		 values_fill = list(alt_aa_freq = 0)
)

## ----doParCoordPlot, fig.width = 4, fig.height = 8, out.width = "50%", results = "hide"----
par(mfcol = c(4,1), mar = c(2, 4, 2, 2) + 0.1)
for(case in sort(unique(wide_ihv$case))){
	thisCase <- t(wide_ihv[wide_ihv$case == case, c(15, 16, 17, 18, 19)])
	matplot(thisCase, type = "b", pch = 19, col = rgb(0, 0, 0, 0.4), lty = 1,
		xaxt = "n", xlab = "",
		main = sub("c", "case ", case),
		ylab = "SNV frequency"
	)
	axis(side = 1, at = 1:5, labels = rownames(thisCase)
	)
}

## ----c18Contrasts-------------------------------------------------------------
intraHostVariants[intraHostVariants$case == "c18" & intraHostVariants$alt_aa_freq > 0.2 ,
                  c(1, 2, 10, 13, 14, 15, 16)]

## ----splitData, results = "hide"----------------------------------------------
strata <- split(intraHostVariants,
                    list(intraHostVariants$organ, intraHostVariants$case),
                    sep = "_", drop = TRUE
)

## ----createGenome, results = "hide"-------------------------------------------
mn908947.3.ihv <- mn908947.3
# Get the nsp12 3' end and assign it to nsp12_end
nsp12_rows <- mn908947.3.ihv$CDS$protein %in% c("nsp12a", "nsp12b")
nsp12_end <- max(mn908947.3.ihv$CDS$end[nsp12_rows])
# Create an entry for nsp12 from nsp12a
mn908947.3.ihv$CDS$protein[mn908947.3.ihv$CDS$protein == "nsp12a"] <- "nsp12"
# Update the 3' end
mn908947.3.ihv$CDS$end[mn908947.3.ihv$CDS$protein == "nsp12"] <- nsp12_end
# dismiss the old nsp12b annotation
mn908947.3.ihv$CDS <- mn908947.3.ihv$CDS[mn908947.3.ihv$CDS$protein != "nsp12b", ]
# tidy up the environment
rm(nsp12_rows, nsp12_end)

## ----computeEntropies, results = "hide"---------------------------------------
profiles <- lapply(strata, function(df) {
			   getEntropySignature(df, position = "POS", ref = "REF", alt = "ALT",
					       genome = mn908947.3.ihv
			   )
})

## ----doHeatMaps, fig.width = 9, fig.height = 7, out.width = "90%", results = "hide"----
heatmap_entropyProfiles(!!!profiles)

## ----doBoxPlots, fig.width = 9, fig.height = 7, out.width = "90%", results = "hide"----
combined_entropy <- lapply(profiles, function(p) p$Entropy) |> do.call(rbind, args = _)
combined_entropy <- cbind(combined_entropy,
                          system = factor(ifelse(grepl("lung", rownames(combined_entropy)),
                                                       "respiratory system", "other"),
                                          levels = c("respiratory system", "other")
                                   )
                    )
# Plot in genomic order
combined_entropy$protein <- factor(
                                   combined_entropy$protein,
                                   levels = mn908947.3.ihv$CDS$protein
)
bp <- ggplot2::ggplot(data = combined_entropy, aes(x = protein, y = entropy)) +
	ggplot2::geom_boxplot(varwidth = T) +
	ggplot2::facet_wrap(~ system, ncol = 1, scales = "free_y") +
	ggplot2::scale_x_discrete(drop = FALSE) +
	ggplot2::theme_bw() +
	ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
bp

