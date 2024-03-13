############################# PARP2_summary ###############################################
data(PARP2_summary)
Breast <- PARP2_summary$Breast
Thyroid <- PARP2_summary$Thyroid
Betah <- list(Breast$beta, Thyroid$beta)
Sigmah <- list(diag(Breast$se), diag(Thyroid$se))
genename <- "PARP2"
snpnames <- Breast$snp
K <- 2
m <- 6
set.seed(123)
RES <- HS(Betah, Sigmah,
  kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
  m = m, K = K, niter = 800, burnin = 400, nthin = 1, nchains = 1,
  a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
)
\dontrun{
  print(RES)
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

  RES <- HS(Betah, Sigmah,
    kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 1, nchains = 1,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
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

  RES <- HS(Betah, Sigmah,
    kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
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


  RES <- HS(Betah, Sigmah,
    kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
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


  RES <- HS(Betah, Sigmah,
    kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
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

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES1)
}
