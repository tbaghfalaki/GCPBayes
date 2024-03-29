% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaryCS.R
\name{summaryCS}
\alias{summaryCS}
\title{Summary function of Continuous Spike}
\usage{
summaryCS(object)
}
\arguments{
\item{object}{a result of a call to CS}
}
\value{
\itemize{
\item Name of Gene: the component from object.
\item Number of SNPs: the component from object.
\item Name of SNPs: the component from object.
\item log10BF: the component from object.
\item lBFDR: the component from object.
\item theta: the component from object.
\item Significance based on CI: the component from object.
}
}
\description{
summaryCS is a generic function used to produce result summaries of the results of CS function.
}
\examples{
############################# Gene DNAJC1 ###############################################
data(DNAJC1)
Breast <- DNAJC1$Breast
Thyroid <- DNAJC1$Thyroid
genename <- "DNAJC1"
snpnames <- Breast$snp
Betah <- list(Breast$beta, Thyroid$beta)
Sigmah <- list(diag(Breast$se^2), diag(Thyroid$se^2))
K <- 2
m <- 14

pvalue <- matrix(0, K, m)
for (k in 1:K) {
  pvalue[k, ] <- 2 * pnorm(-abs(Betah[[k]] / sqrt(diag(Sigmah[[k]]))))
}

zinit <- rep(0, K)
for (j in 1:K) {
  index <- 1:m
  PVALUE <- p.adjust(pvalue[j, ])
  SIGNALS <- index[PVALUE < 0.05]
  modelf1 <- rep(0, m)
  modelf1[SIGNALS] <- 1
  if (max(modelf1) == 1) (zinit[j] <- 1)
}


RES <- CS(Betah, Sigmah,
  kappa0 = 0.5, tau20 = 1, zeta0 = zinit,
  m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1, a1 = 0.1, a2 = 0.1,
  c1 = 0.1, c2 = 0.1, sigma2 = 10^-3, snpnames = snpnames, genename = genename
)
summaryCS(RES)
}
\references{
Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. Statistics in Medicine, 40(6), 1498-1518.
}
\author{
Taban Baghfalaki.
}
