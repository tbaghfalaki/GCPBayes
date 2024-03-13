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

set.seed(123)
RES <- DS(Betah, Sigmah,
  kappa0 = 0.5, sigma20 = 1,
  m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1,
  a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, snpnames, genename
)
summaryDS(RES)
