#' SNVs from multiple organs of individuals with postmortem SARS-CoV-2 detection
#' 
#' SNVs were inferred from Illumina (2 x 150) sequences. Reads were aligned to
#' the Wuhan-Hu-1 reference genome (MN908947.3) with \code{BWA}, followed by
#' variant calling using \code{LoFreq} with source quality correction. The VCF
#' file was post-processed with \code{BCFtools} to prevent reference and strand
#' biases and remove positions supported by fewer than 1,000 reads.
#'
#' @docType data
#' @format A data frame with 213 rows and 16 columns.
#'
#' @seealso 
#' See \code{vignette("intraHostVariants", package = "MetaEntropy")} for a 
#' step-by-step tutorial on analyzing this dataset.
#'
#' @source
#' Manrique, J. M., Maffia-Bizzozero, S., Delpino, M. V., Quarleri, J., &
#' Jones, L. R. (2024). Multi-Organ Spread and Intra-Host Diversity of
#' SARS-CoV-2 Support Viral Persistence, Adaptation, and a Mechanism That
#' Increases Evolvability. \emph{Journal of Medical Virology}, 96(12), e70107. 
#' \doi{10.1002/jmv.70107}
#' 
"intraHostVariants"
