#############################Gene DNAJC1 ###############################################
data(DNAJC1)
Breast <- DNAJC1$Breast
Thyroid <- DNAJC1$Thyroid
genename <- "DNAJC1"
snpnames <- Breast$snp
Betah <- list(Breast$beta, Thyroid$beta)
Sigmah <- list(diag(Breast$se^2), diag(Thyroid$se^2))
K <- 2
m <- 14


RES1 <- DS(Betah, Sigmah,
           kappa0 = 0.5, sigma20 = 1,
           m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1,
           a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, snpnames, genename
)


MCMCplot(Result = RES1, k = 2, nchains = 1, whichsnps = sample(snpnames, 7),
                     betatype = "l",
                     acftype = "correlation",
                     dencol = "white", denlty = 1, denbg = "white")
###################Simulated summary level data with K=5 ###############################
\dontrun{
data(Simulated_summary)
genename <- Simulated_summary$genename
snpnames <- Simulated_summary$snpnames
Betah <- Simulated_summary$simBeta
Sigmah <- Simulated_summary$simSIGMA
K <- 5
m <- 10

RES1 <- DS(Betah, Sigmah,
 kappa0 = c(0.2, 0.5), sigma20 = c(1, 2),
 m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
 a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, snpnames, genename)

MCMCplot(Result = RES1, k = 3, nchains = 2, whichsnps = sample(snpnames, 3),
         betatype = "l",
         acftype = "partial",
         dencol = "blue", denlty = 1, denbg = "black")

RES1 <- DS(Betah, Sigmah,
 kappa0 = c(0.2, 0.5, 0.6), sigma20 = c(1, 2, 1.5),
 m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 3,
 a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, snpnames, genename)

MCMCplot(Result = RES1, k = 3, nchains = 3, whichsnps = sample(snpnames, 5),
         betatype = "l",
         acftype = "partial",
         dencol = "white", denlty = 1, denbg = "white")
#############################Gene DNAJC1 ###############################################
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
c1 = 0.1, c2 = 0.1, sigma2 = 10^-3, snpnames = snpnames, genename = genename)

MCMCplot(Result = RES1, k = 1, nchains = 1, whichsnps = sample(snpnames, 7),
         betatype = "l",
         acftype = "correlation",
         dencol = "white", denlty = 1, denbg = "white")
###################################Gene PARP2 ##########################################
library(BhGLM)
data(PARP2)
Breast <- PARP2$Breast
Thyroid <- PARP2$Thyroid
genename <- "PARP2"
snpnames <- c("rs3093872", "rs3093921", "rs1713411", "rs3093926", "rs3093930", "rs878156")


Fit1 <- BhGLM::bglm(y1~ ., family=binomial(link="logit"),data=Breast)
Betah1 <-  Fit1$coefficients[-1]
Sigmah1 <- cov(coef(arm::sim(Fit1)))[-1,-1]

Fit2 <- BhGLM::bglm(y2~ ., family=binomial(link="logit"),data=Thyroid)
Betah2 <-  Fit2$coefficients[-1]
Sigmah2 <- cov(coef(arm::sim(Fit2)))[-1,-1]

Betah <- list(Betah1,Betah2)
Sigmah <- list(Sigmah1,Sigmah2)
K <- 2
m <- 6


RES1 <- DS(Betah, Sigmah, kappa0=c(0.2,0.5), sigma20=c(1,2),
          m=m, K=K, niter=1000, burnin=500, nthin=1, nchains=2,
          a1=0.1, a2=0.1, d1=0.1, d2=0.1, snpnames, genename)

MCMCplot(Result=RES1, k=1, nchains=2, whichsnps=snpnames,
         betatype = "l",
         acftype = "correlation",
         dencol = "red", denlty = 1, denbg = "white")


RES1 <- DS(Betah, Sigmah, kappa0=c(0.2,0.5), sigma20=c(1,2),
          m=m, K=K, niter=2000, burnin=1000, nthin=2, nchains=2,
          a1=0.1, a2=0.1, d1=0.1, d2=0.1, snpnames, genename)

MCMCplot(Result=RES1, k=1, nchains=2, whichsnps=snpnames,
         betatype = "l",
         acftype = "correlation",
         dencol = "white", denlty = 1, denbg = "white")
}

