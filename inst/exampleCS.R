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

set.seed(123)
RES <- CS(Betah, Sigmah,
  kappa0 = 0.5, tau20 = 1, zeta0 = zinit,
  m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1, a1 = 0.1, a2 = 0.1,
  c1 = 0.1, c2 = 0.1, sigma2 = 10^-3, snpnames = snpnames, genename = genename
)
\dontrun{
  print(RES)


  RES1 <- CS(Betah, Sigmah,
    kappa0 = c(0.2, 0.5), tau20 = c(1, 2), zeta0 = zinit,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, c1 = 0.1, c2 = 0.1, sigma2 = 10^-3, snpnames, genename
  )
  print(RES1)
  ################### Simulated summary level data with K=5 ###############################
  data(Simulated_summary)
  genename <- Simulated_summary$genename
  snpnames <- Simulated_summary$snpnames
  Betah <- Simulated_summary$simBeta
  Sigmah <- Simulated_summary$simSIGMA
  K <- 5
  m <- 10

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
  print(RES)


  RES1 <- CS(Betah, Sigmah,
    kappa0 = c(0.2, 0.5), tau20 = c(1, 2), zeta0 = zinit,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, c1 = 0.1, c2 = 0.1, sigma2 = 10^-3, snpnames, genename
  )
  print(RES1)
  ################################### Gene PARP2 ##########################################
  library(BhGLM)
  data(PARP2)
  Breast <- PARP2$Breast
  Thyroid <- PARP2$Thyroid
  genename <- "PARP2"
  snpnames <- c("rs3093872", "rs3093921", "rs1713411", "rs3093926", "rs3093930", "rs878156")


  Fit1 <- BhGLM::bglm(y1 ~ ., family = binomial(link = "logit"), data = Breast)
  Betah1 <- Fit1$coefficients[-1]
  Sigmah1 <- cov(coef(arm::sim(Fit1)))[-1, -1]

  Fit2 <- BhGLM::bglm(y2 ~ ., family = binomial(link = "logit"), data = Thyroid)
  Betah2 <- Fit2$coefficients[-1]
  Sigmah2 <- cov(coef(arm::sim(Fit2)))[-1, -1]

  Betah <- list(Betah1, Betah2)
  Sigmah <- list(Sigmah1, Sigmah2)
  K <- 2
  m <- 6

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
  print(RES)

  RES1 <- CS(Betah, Sigmah,
    kappa0 = c(0.2, 0.5), tau20 = c(1, 2), zeta0 = zinit,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, c1 = 0.1, c2 = 0.1, sigma2 = 10^-3, snpnames, genename
  )
  print(RES1)
  ########### Simulated individual level data with K=3 and continuous phynotype ###########
  library(BhGLM)
  data(Simulated_individual)
  Study1 <- Simulated_individual$Study1
  Study2 <- Simulated_individual$Study2
  Study3 <- Simulated_individual$Study3
  K <- 3
  m <- 30
  genename <- "Simulated"
  snpnames <- sprintf("SNP%s", seq(1:m))


  Fit1 <- BhGLM::bglm(Y1 ~ ., family = gaussian, data = data.frame(Study1))
  Betah1 <- Fit1$coefficients[-1]
  Sigmah1 <- cov(coef(arm::sim(Fit1)))[-1, -1]

  Fit2 <- BhGLM::bglm(Y2 ~ ., family = gaussian, data = data.frame(Study2))
  Betah2 <- Fit2$coefficients[-1]
  Sigmah2 <- cov(coef(arm::sim(Fit2)))[-1, -1]

  Fit3 <- BhGLM::bglm(Y3 ~ ., family = gaussian, data = data.frame(Study3))
  Betah3 <- Fit3$coefficients[-1]
  Sigmah3 <- cov(coef(arm::sim(Fit3)))[-1, -1]

  Betah <- list(Betah1, Betah2, Betah3)
  Sigmah <- list(Sigmah1, Sigmah2, Sigmah3)


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
  print(RES)

  RES1 <- CS(Betah, Sigmah,
    kappa0 = c(0.2, 0.5), tau20 = c(1, 2), zeta0 = zinit,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, c1 = 0.1, c2 = 0.1, sigma2 = 10^-3, snpnames, genename
  )
  print(RES1)

  ########### Simulated individual level data with K=2 and gene expression data ###########
  library(BhGLM)
  data(Simulated_individual_survival)
  Study1 <- Simulated_individual_survival$Study1
  Study2 <- Simulated_individual_survival$Study2
  K <- 2
  m <- 10
  genename <- "Simulated"
  snpnames <- sprintf("G%s", seq(1:m))


  Fit1 <- BhGLM::bcoxph(Study1$T ~ Study1$X)
  Betah1 <- Fit1$coefficients
  Sigmah1 <- Fit1$var


  Fit2 <- BhGLM::bcoxph(Study2$T ~ Study2$X)
  Betah2 <- Fit2$coefficients
  Sigmah2 <- Fit2$var

  Betah <- list(Betah1, Betah2)
  Sigmah <- list(Sigmah1, Sigmah2)

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


  RES1 <- CS(Betah, Sigmah,
    kappa0 = c(0.2, 0.5), tau20 = c(1, 2), zeta0 = zinit,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, c1 = 0.1, c2 = 0.1, sigma2 = 10^-3, snpnames, genename
  )
  print(RES1)
}
