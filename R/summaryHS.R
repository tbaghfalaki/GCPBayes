#' Summary function of Hierarchical Spike
#'
#'
#' @description
#' The summaryHS function serves as a versatile tool for generating result summaries derived from the output of the HS function.
#'
#'
#'
#'
#' @param object a result of a call to HS
#'
#' @return
#' - Name of Gene: the component from object.
#' - Number of SNPs: the component from object.
#' - Name of SNPs: the component from object.
#' - Pleiotropic effect based on CI: the component from object.
#' - Pleiotropic effect based on median thresholding: the component from object.
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. Statistics in Medicine, 40(6), 1498-1518.
#'
#' @example inst/examplesummaryHS.R
#'
#' @md
#' @export



summaryHS <- function(object) {
  list(
    "Name of Gene" = object$Criteria$`Name of Gene`,
    "Number of SNPs" = length(object$Criteria$`Name of SNPs`),
    "Name of SNPs" = object$Criteria$`Name of SNPs`,
    "Pleiotropic effect based on CI" = object$Indicator$`Significant studies and Variable pleiotropic effect based on CI`$`Pleiotropic effect`,
    "Pleiotropic effect based on median thresholding" = object$Indicator$`Significant studies and Variable pleiotropic effect based on median`$`Pleiotropic effect`
  )
}
