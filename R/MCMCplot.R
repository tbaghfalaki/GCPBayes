#'  MCMC plot
#'
#'
#' @description
#'  Trace plot, density plot and ACF plot for the output of CS/DS/HS. The plot is able to draw at most ten SNPs.
#'
#'
#' @details
#'  Trace plot, density plot and ACF plot for the output of CS/DS/HS for checking convergence of MCMC chains.
#'
#'
#' @param Result All the generated results by CS/DS/HS function.
#' @param k The number of study for drawing plots, k=1,2,...,K.
#' @param nchains Number of Markov chains run in Result.
#' @param whichsnps The name of SNPs.
#' @param betatype The type of plot desired. The following values are possible: "p" for points, "l" for lines, "b" for both points and lines, "c" for empty points joined by lines, "o" for overplotted points and lines, "s" and "S" for stair steps and "h" for histogram-like vertical lines. Finally, "n" does not produce any points or lines.
#' @param acftype String giving the type of ACF to be computed. Allowed values are "correlation" (the default), "covariance" or "partial". Will be partially matched.
#' @param dencol The color for filling the density plot.
#' @param denlty The line type to be used in the density plot.
#' @param denbg The color to be used for the background of the density plot.
#'
#' @importFrom graphics lines par
#' @importFrom stats ACF density
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. Statistics in Medicine, 40(6), 1498-1518.
#'
#' @example inst/exampleMCMCplot.R
#'
#' @md
#' @export

MCMCplot <- function(Result = Result, k = k, nchains = nchains, whichsnps = whichsnps,
                     betatype = "l",
                     acftype = "correlation",
                     dencol = "white", denlty = 1, denbg = "white") {
  jjj <- k
  RES1 <- Result
  whichones <- whichsnps
  snpnames <- Result$Criteria$`Name of SNPs`
  if (length(whichones) <= 10) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(3, length(whichones)))

    A <- RES1$MCMCChain[[1]]$Beta[[jjj]]
    colnames(A) <- snpnames
    if (nchains == 2) {
      B <- RES1$MCMCChain[[2]]$Beta[[jjj]]
      colnames(B) <- snpnames
    }
    if (nchains > 2) {
      B <- RES1$MCMCChain[[2]]$Beta[[jjj]]
      colnames(B) <- snpnames
      C <- RES1$MCMCChain[[3]]$Beta[[jjj]]
      colnames(C) <- snpnames
    }
    for (k in 1:length(whichones)) {
      P0 <- plot(A[, whichones][, k],
        type = betatype,
        main = whichones[k], ylab = paste("beta", k)
      )
      # print(P0)
      if (nchains == 2) {
        L0 <- lines(B[, whichones][, k], col = 2)
        # print(L0)
      }
      if (nchains > 2) {
        L1 <- lines(B[, whichones][, k], col = 2)
        # print(L1)
        L2 <- lines(C[, whichones][, k], col = 3)
        # print(L2)
      }
    }

    for (k in 1:length(whichones)) {
      den <- density(A[, whichones][, k])
      # Create basic plot
      plot(den, main = "")
      # Change the plot region color
      rect(par("usr")[1], par("usr")[3],
        par("usr")[2], par("usr")[4],
        col = denbg
      )
      # Add a new plot
      par(new = TRUE)
      plot(den, main = whichones[k])
      polygon(den, col = dencol, lty = denlty)
      # print(P1)
    }
    for (k in 1:length(whichones)) {
      ACF <- acf(A[, whichones][, k],
        main = whichones[k],
        type = acftype
      )
      # print(ACF)
    }
  } else {
    # print("Please consider at most 10 variables/SNPs.")
  }
}
