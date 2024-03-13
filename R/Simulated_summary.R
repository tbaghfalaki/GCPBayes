#' Simulated summary statistics for K=5 traits
#'
#' A list containing the summary statistics including regression coefficients and covariance matrices for K=5 studies.
#' @name Simulated_summary
#' @format A list which contains the name of gene, the name of the SNPs, the vectors and the matrices for the summary statistics of each study.
#' \describe{
#'   \item{genename}{The name of gene}
#'   \item{snpnames}{The name of the SNPs}
#'   \item{simBeta}{The regression coeffiecnts of the studies}
#'   \item{simSIGMA}{The covariance matrices for the studies}
#' }
#' @references
#' \enumerate{'
#' \item
#' Baghfalaki, T., Sugier, Y. Asgari, P. E., Truong, & Liquet, B. (2021). GCPBayes:  An R Package for Studying Cross-Phenotype Genetic Associations with Group-Level Bayesian Meta-Analysis. \emph{Submitted}.
#' }
#' @seealso \code{\link{GCPBayes}}
"Simulated_summary"
