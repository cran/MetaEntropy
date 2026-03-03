## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----loadMetaEntropy, results = "hide"----------------------------------------
library(MetaEntropy)

## ----seeLinked----------------------------------------------------------------
wWater[105:108,-1]

## ----seewMultiple-------------------------------------------------------------
wWater[50:51,-1]

## ----loadPackages, results = "hide"-------------------------------------------
lapply(c("ggplot2", "patchwork"), library, character.only = TRUE)

## ----datasets-----------------------------------------------------------------
firstWave <- wWater[ wWater$wave == "first", ]
thirdWave <- wWater[ wWater$wave == "third", ]

## ----filter-------------------------------------------------------------------
firstWave <- firstWave[ firstWave$alt_aa_freq <= 0.97, ]
thirdWave <- thirdWave[ thirdWave$alt_aa_freq <= 0.97, ]

## ----signatures---------------------------------------------------------------
ancestral <- getEntropySignature(firstWave)
omicron <-  getEntropySignature(thirdWave)
# Compare signatures
anc_plot <- plot(ancestral) + ggtitle("Ancestral")
omi_plot <- plot(omicron) + ggtitle("Omicron")
anc_plot / omi_plot

## ----omicronEntroScan, fig.width = 7------------------------------------------
plot(omicron, chartType = "entroScan")

## -----------------------------------------------------------------------------
omicron$Entropy$position[ omicron$Entropy$entropy > 0.3 ]

## -----------------------------------------------------------------------------
showMutations(omicron, c(22882, 22898, 22917, 23013, 23040, 23048, 23055, 23063))

## ----eval = FALSE-------------------------------------------------------------
# # Search PUBMED for more studies on these mutations:
# library(rentrez)
# #
# # create search phrase
# mutations <- showMutations(omicron, c(22882, 22898, 22917, 23013, 23040, 23048, 23055, 23063))
# searchPhrase <- paste0(gsub("S:", "", mutations$phenotype, fixed = TRUE), "[TIAB]")
# searchPhrase <- paste0("SARS-CoV-2[TITLE] AND spike [TIAB] AND (", paste0(searchPhrase, collapse = " OR "), ")")
# #
# # search
# mySearch <- entrez_search(db = "pubmed", term = searchPhrase, retmax = 1000)
# #
# # Number of studies
# mySearch$count
# [1] 560
# #

## -----------------------------------------------------------------------------
assessHotSpot(omicron, c(22517, 23186))

## -----------------------------------------------------------------------------
assessHotSpot(ancestral, c(22517, 23186))

## -----------------------------------------------------------------------------
summary(ancestral)

## -----------------------------------------------------------------------------
summary(omicron)

