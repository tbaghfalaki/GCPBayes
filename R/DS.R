#' Dirac Spike
#'
#'
#' @description
#' Utilize a Gibbs sampler to conduct analysis on a multivariate Bayesian sparse group selection model, employing a Dirac spike prior to identify pleiotropic effects across traits. This function is tailored for summary statistics comprising estimated regression coefficients and their corresponding covariance matrices.
#'
#'
#' @details
#' Let betah_k, k=1,...,K be a m-dimensional vector of the regression coefficients for the kth study and Sigmah_k  be its estimated covariance matrix. The hierarchical set-up of  DS prior, by considering summary statistics (betah_k and Sigmah_k, k=1,...,K) as the input of the method, is given by:
#'
#' betah _k ~ (1 - kappa) delta_0(betah_k) + kappa N_m(0,sigma2 I_m ),
#'
#' kappa ~ Beta(a_1,a_2),
#'
#' sigma2 ~ inverseGamma (d_1,d_2).
#'
#' where delta_0(betah_k) denotes a point mass at 0, such that delta_0(betah_k)=1 if beta_k=0 and  delta_0(betah_k)=0 if  at least one of the $m$ components of beta_k is non-zero.
#'
#'
#'
#' @param Betah A list containing m-dimensional vectors of the regression coefficients for K studies.
#' @param Sigmah A list containing the positive definite covariance matrices (m*m-dimensional) which is the estimated covariance matrices of K studies.
#' @param kappa0 Initial value for kappa (its dimension is equal to nchains).
#' @param sigma20 Initial value for sigma2 (its dimension is equal to nchains).
#' @param m Number of variables in the group.
#' @param K  Number of traits.
#' @param niter Number of iterations for the Gibbs sampler.
#' @param burnin Number of burn-in iterations.
#' @param nthin The lag of the iterations used for the posterior analysis is defined (or thinning rate).
#' @param nchains Number of Markov chains, when nchains>1, the function calculates the Gelman-Rubin convergence statistic, as modified by Brooks and Gelman (1998).
#' @param a1,a2 Hyperparameters of kappa. Default is a1=0.1 and a2=0.1.
#' @param d1,d2 Hyperparameters of sigma2. Default is d1=0.1 and d2=0.1.
#' @param snpnames Names of variables for the group.
#' @param genename Name of group.
#'
#' @importFrom stats median quantile rbeta rbinom sd
#'
#' @return
#' - mcmcchain: The list of simulation output for all parameters.
#' - Summary: Summary statistics for regression coefficients in each study.
#' - Criteria: genename, snpnames, PPA, log10BF, lBFDR, theta.
#' - Indicator: A table containing m rows of binary indicators for each study, the number of studies with nonzero signal and having pleiotropic effect by credible interval (CI). The first K columns show nonzero signals, K+1 th column includes the number of studies with nonzero signal and the last column shows an indicator for having pleiotropic effect of each SNP.
#'
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' \enumerate{
#' \item
#' Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. \emph{Statistics in Medicine}, \strong{40}(6), 1498-1518.
#' }
#' @example inst/exampleDS.R
#'
#' @md
#' @export


DS <- function(Betah, Sigmah, kappa0, sigma20, m, K, niter = 1000, burnin = 500, nthin = 2, nchains = 2, a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, snpnames, genename) {
  RES1 <- list()
  Result <- list()
  for (j in 1:nchains) {
    RES1[[j]] <- DS0(Betah, Sigmah,
      kappa0 = kappa0[j], sigma20 = sigma20[j], m = m,
      K = K, niter = niter, burnin = burnin, nthin = nthin, a1 = a1, a2 = a2, d1 = d1, d2 = d2,
      snpnames = snpnames, genename = genename
    )

    Result[[j]] <- RES1[[j]]$mcmcchain
  }

  Summary <- list()
  if (m > 1) {
    if (nchains > 1) {
      ts <- sample(1:nchains, 2)
      for (k in 1:K) {
        AA <- list(RES1[[ts[1]]]$mcmcchain$Beta[[k]], RES1[[ts[2]]]$mcmcchain$Beta[[k]]) # Study k
        Summary$Beta[[k]] <- data.frame(
          snpnames, apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, mean), apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, sd),
          t(apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, function(x) quantile(x, c(.025, 0.5, .975)))),
          wiqid::simpleRhat(postpack::post_convert(AA))
        )
        colnames(Summary$Beta[[k]]) <- cbind("Name of SNP", "Mean", "SD", "val2.5pc", "Median", "val97.5pc", "BGR")
      }
    }
    if (nchains == 1) {
      for (k in 1:K) {
        Summary$Beta[[k]] <- data.frame(
          snpnames, apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, mean), apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, sd),
          t(apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, function(x) quantile(x, c(.025, 0.5, .975))))
        )
        colnames(Summary$Beta[[k]]) <- cbind("Name of SNP", "Mean", "SD", "val2.5pc", "Median", "val97.5pc")
      }
    }
  }
  if (m == 1) {
    if (nchains > 1) {
      ts <- sample(1:nchains, 2)
      for (k in 1:K) {
        AA <- list(matrix(RES1[[ts[1]]]$mcmcchain$Beta[[k]]), matrix(RES1[[ts[2]]]$mcmcchain$Beta[[k]])) # Study k
        Summary$Beta[[k]] <- data.frame(
          snpnames, mean(RES1[[1]]$mcmcchain$Beta[[k]]), sd(RES1[[1]]$mcmcchain$Beta[[k]]),
          t(quantile(RES1[[1]]$mcmcchain$Beta[[k]], c(.025, 0.5, .975))),
          wiqid::simpleRhat(postpack::post_convert(AA))
        )
        colnames(Summary$Beta[[k]]) <- cbind("Name of SNP", "Mean", "SD", "val2.5pc", "Median", "val97.5pc", "BGR")
      }
    }
    if (nchains == 1) {
      for (k in 1:K) {
        Summary$Beta[[k]] <- data.frame(
          snpnames, mean(RES1[[1]]$mcmcchain$Beta[[k]]), sd(RES1[[1]]$mcmcchain$Beta[[k]]),
          t(quantile(RES1[[1]]$mcmcchain$Beta[[k]], c(.025, 0.5, .975)))
        )
        colnames(Summary$Beta[[k]]) <- cbind("Name of SNP", "Mean", "SD", "val2.5pc", "Median", "val97.5pc")
      }
    }
  }



  RES1new <- list(MCMCChain = Result, Summary = Summary, Criteria = RES1[[1]]$Criteria, Indicator = RES1[[1]]$Indicator)


  return(RES1new)
}


DS0 <- function(Betah, Sigmah, kappa0, sigma20, m, K, niter = 100, burnin, nthin = 1, a1 = 0.1, a2 = 1, d1 = 0.1, d2 = 1, snpnames, genename) {
  # memory.limit(99999999)

  Beta <- list()
  for (k in 1:K) {
    Beta$Beta[[k]] <- matrix(0, niter, m)
  }

  for (k in 1:K) {
    Beta$Beta[[k]][1, ] <- Betah[[k]]
  }


  kappa <- rep(0, niter)
  kappa[1] <- kappa0
  sigma2 <- rep(0, niter)
  sigma2[1] <- sigma20

  psi1 <- psi2 <- rep(1, niter)

  MeanBeta1 <- MeanBeta2 <- MedianBeta1 <- MedianBeta2 <- c()
  Geneplotci <- Geneplotmed <- c()
  PROB1j <- PROB2j <- c()
  PROB1 <- PROB2 <- c()
  Psi <- matrix(0, K, niter)
  for (r in 2:niter) {
    ##################### Beta_k ####
    for (k in 1:K) {
      Sigmahk <- Sigmah[[k]] # Sigmah[((k-1)*m+1):(k*m),((k-1)*m+1):(k*m)]
      Omega1 <- arma_inv(Sigmahk) + (1 / sigma2[r - 1]) * diag(m)
      Mean1 <- arma_inv(Omega1) %*% arma_inv(Sigmahk) %*% matrix(Betah[[k]], m, 1)
      kappat1 <- 1 / (1 + (kappa[r - 1] / (1 - kappa[r - 1]) *
        exp(.5 * t(Mean1) %*% Omega1 %*% Mean1 - .5*arma_log_det(Omega1) - m / 2 * log(sigma2[r - 1]))))

      Tab <- rbinom(1, 1, as.numeric(kappat1))
      if (Tab == 0) (Beta$Beta[[k]][r, ] <- mvtnorm::rmvnorm(1, mean = Mean1, sigma = arma_inv(Omega1)))
      if (Tab == 1) (Beta$Beta[[k]][r, ] <- rep(0, m))
      Psi[k, r] <- 1 - Tab
    }

    ##################### kappa ####
    K0 <- 0
    for (k in 1:K) {
      if (max(Beta$Beta[[k]][r, ]) == 0) (K0 <- K0 + 1)
    }
    kappa[r] <- rbeta(1, K - K0 + a1, K0 + a2)

    ##################### sigma2 ####
    tBB <- c()
    for (k in 1:K) {
      tBB[k] <- t(Beta$Beta[[k]][r, ]) %*% Beta$Beta[[k]][r, ]
    }
    sigma2[r] <- invgamma::rinvgamma(1,
      shape = m * (K - K0) / 2 + d1,
      rate = d2 + sum(tBB) / 2
    )
  }

  # $$$$$$$$$$$$$$$$$$$$$$$$$ New Posterior samples After Burn-in $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  niter1 <- burnin + 1
  index <- niter1:niter
  indexn <- index[(index %% nthin) == 0]


  sigma2 <- sigma2[indexn]
  kappa <- kappa[indexn]
  Psi <- Psi[, indexn]
  PPA <- list()
  for (k in 1:K) {
    Beta$Beta[[k]] <- Beta$Beta[[k]][indexn, ]
    PPA[[k]] <- mean(Psi[k, ])
  }
  Beta0 <- Beta$Beta
  PGENE <- PGENE1 <- matrix(0, m, K)
  B1CI <- B2CI <- c()
  if (m > 1) {
    for (k in 1:K) {
      for (int in 1:m) {
        B1CI <- quantile(Beta$Beta[[k]][, int], c(0.025, 0.975))
        if ((0 < B1CI[1]) | (0 > B1CI[2])) (PGENE[int, k] <- 1)
      }
    }
  } else {
    for (k in 1:K) {
      B1CI <- quantile(Beta$Beta[[k]], c(0.025, 0.975))
      if ((0 < B1CI[1]) | (0 > B1CI[2])) (PGENE[1, k] <- 1)
    }
  }


  Geneplotci <- apply(PGENE, 1, sum)


  pe <- rep("No", m)
  total <- apply(PGENE, 1, sum)
  pe[total > 1] <- "Yes"
  Geneplotci <- data.frame(PGENE, total, pe)
  rownames(Geneplotci) <- snpnames
  colnames(Geneplotci) <- c(paste("Study", 1:K), "Total", "Pleiotropic effect")


  #####################################
  if (m > 1) {
    for (k in 1:K) {
      for (int in 1:m) {
        B1med <- quantile(Beta$Beta[[k]][, int], c(0.5))
        if (B1med != 0) (PGENE1[int, k] <- 1)
      }
    }
  } else {
    for (k in 1:K) {
      B1med <- quantile(Beta$Beta[[k]], c(0.5))
      if (B1med != 0) (PGENE1[1, k] <- 1)
    }
  }


  Geneplotmed <- apply(PGENE1, 1, sum)
  pemed <- rep("No", m)
  total1 <- apply(PGENE1, 1, sum)
  pemed[total1 > 1] <- "Yes"
  Geneplotmedian <- data.frame(PGENE1, total1, pemed)
  rownames(Geneplotmedian) <- snpnames
  colnames(Geneplotmedian) <- c(paste("Study", 1:K), "Total", "Pleiotropic effect")

  # $$$$$$$$$$$$$$$$$$$$$$$$ locFDR, theta, BF $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  pzeta <- matrix(0, length(indexn), K)

  pz <- theta <- rep(0, length(indexn))
  for (r in 1:(length(indexn))) {
    for (k in 1:K) {
      Sigmahk <- Sigmah[[k]] # [((k-1)*m+1):(k*m),((k-1)*m+1):(k*m)]
      Omega1 <- arma_inv(Sigmahk) + (1 / sigma2[r]) * diag(m)
      Mean1 <- arma_inv(Omega1) %*% arma_inv(Sigmahk) %*% matrix(Betah[[k]], m, 1)

      pzeta[r, k] <- 1 / (1 + (kappa[r] / (1 - kappa[r]) *
        exp(.5 * t(Mean1) %*% Omega1 %*% Mean1 - .5 *arma_log_det(Omega1) - m / 2 * log(sigma2[r]))))
    }

    pz[r] <- prod(pzeta[r, 1:K])

    if (K != 2) {
      except1 <- rep(0, K)
      for (j in 2:(K - 1)) {
        except1[j] <- prod((1 - pzeta[r, j]) * (pzeta[r, c(1:(j - 1), (j + 1):K)]))
      }
      theta[r] <- 1 - prod(pzeta[r, 1:K]) -
        prod((1 - pzeta[r, 1]) * (pzeta[r, c(2:K)])) - sum(except1[-1]) -
        prod(pzeta[r, 1:(K - 1)] * (1 - pzeta[r, K]))
    } else {
      theta[r] <- (1 - pzeta[r, 1]) * (1 - pzeta[r, 2])
    }
  }


  locFDR <- mean(pz[is.nan(pz) == FALSE])
  locFDR
  ResFDR <- locFDR

  BF <- ((1 - locFDR) / (locFDR)) * (a2^K / ((a1 + a2)^K - a2^K))
  BF
  if (locFDR == 0) (BF <- ((1 - locFDR) / (locFDR + 10^-30)) * (a2^K / ((a1 + a2)^K - a2^K)))
  Pleop <- mean(theta[is.nan(theta) == FALSE])
  ResPleop <- Pleop

  ResBF <- log10(BF)
  mcmcchain <- list(kappa = kappa, sigma2 = sigma2, Beta = Beta0)
  Indicator <- list(
    "Significant studies and Pleiotropic effect based on CI" = Geneplotci,
    "Significant studies and Pleiotropic effect based on median thresholding" = Geneplotmedian
  )
  Reslast <- list(
    "Name of Gene" = genename, "Name of SNPs" = snpnames, "PPA" = PPA,
    log10BF = ResBF, lBFDR = ResFDR, theta = ResPleop
  )
  list(mcmcchain = mcmcchain, Criteria = Reslast, Indicator = Indicator)
}
