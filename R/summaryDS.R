#' Summary function of Dirac Spike
#'
#'
#' @description
#' The summaryDS function serves as a versatile tool for generating result summaries derived from the output of the DS function.
#'
#'
#'
#'
#' @param object a result of a call to DS
#'
#' @return
#' - Name of Gene: the component from object.
#' - Number of SNPs: the component from object.
#' - Name of SNPs: the component from object.
#' - log10BF: the component from object.
#' - lBFDR: the component from object.
#' - theta: the component from object.
#' - Significance based on CI: the component from object.
#' - Significance based on median thresholding: the component from object.
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. Statistics in Medicine, 40(6), 1498-1518.
#'
#' @example inst/examplesummaryDS.R
#'
#' @md
#' @export



summaryDS <- function(object) {
  list(
    "Name of Gene" = object$Criteria$`Name of Gene`,
    "Number of SNPs" = length(object$Criteria$`Name of SNPs`),
    "Name of SNPs" = object$Criteria$`Name of SNPs`,
    "log10BF" = object$Criteria$log10BF,
    "lBFDR" = object$Criteria$lBFDR,
    "theta" = object$Criteria$theta,
    "Significance based on CI" = object$Indicator$`Significant studies and Pleiotropic effect based on CI`,
    "Significance based on median thresholding" = object$Indicator$`Significant studies and Pleiotropic effect based on median thresholding`
  )
}
