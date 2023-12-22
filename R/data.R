#' Normalized expression data from TCGA BRCA samples
#'
#'
#' @format ## `TCGA.rsem`
#' A data frame with 48 genes (rows) and 674 samples (columns):
#' \describe{
#'   \item{genes}{Genes useful to classify samples with centroids.}
#'   \item{samples}{Breast cancer samples from the TCGA consortium filtered
#'   for tumor stage II and III}
#'   \item{expression}{The values are RSEM normalized}
#' }
#' @source <https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz>
"TCGA.rsem"
