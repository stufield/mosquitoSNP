#' Genetic analysis of mosquito SNP data
#'
#' A suite of functions & data developed by Stu Field and
#' Martin Donnelly during the Fall 2010 for the analysis of genetic data
#' of insecticide resistant mosquito populations in Africa.
#'
#' Some functionality may depend the \pkg{limma} package,
#' which is suggested but not required.
#' Limma must be downloaded from BioConductor:\cr
#'
#' \verb{
#' if ( !requireNamespace("BiocManager", quietly = TRUE) ) {
#'   install.packages("BiocManager")
#' }
#' BiocManager::install("limma")
#' }
#'
#' @section Borrowed from \pkg{pegas}:
#' Some functionality in this package is derived from the
#' \pkg{pegas} package by Emmanuel Paradis
#' (\url{http://ape-package.ird.fr/pegas.html}) in e.g.
#' calculating heterozygosities. To simplify installation and buffer
#' against version changes in the dependency, this functionality was included
#' as non-exported code and included in a file `pegas.R`.
#'
#' @name mosquitoSNP-package
#' @aliases mosquitoSNP-package mosquito
#' @docType package
#' @seealso Package Data:
#'
#' \code{\link{Ago2}} \cr
#' \code{\link{ambiguityCodons}} \cr
#' \code{\link{codons}} \cr
#' \code{\link{Ghana_MvS_SNPdata}} \cr
#' \code{\link{GhanaKDR}} \cr
#' \code{\link{HarrCraigData}} \cr
#' \code{\link{KenyaSNP_EHHdata}} \cr
#' \code{\link{microsatdata}} \cr
#' \code{\link{Tororo}} (RG.corrected) \cr
#'
#' @references
#' Martin J. Donnelly \cr
#' Vector Group \cr
#' Liverpool School of Tropical Medicine \cr
#' Liverpool L35QA \cr
#' England, UK \cr
#'
#' Stu Field \cr
#' Department of Biology \cr
#' Colorado State University \cr
#' Fort Collins, CO  80523-1878 \cr
#'
#' Bill Black IV \cr
#' Department of Microbiology, Immunology, and Pathology \cr
#' Arthropod-borne and Infectious Diseases Laboratory \cr
#' Colorado State University \cr
#' Fort Collins, CO  80523 \cr
#'
#' @keywords package
#' @examples
#' \dontrun{
#' ambiguityCodons
#' codons
#' head(HarrCraigData)
#' }
NULL
