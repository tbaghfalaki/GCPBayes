#' Hierarchical Spike
#'
#'
#' @description
#' Run a Gibbs sampler for a multivariate Bayesian sparse group selection model with hierarchical spike prior for detecting pleiotropic effects on two traits. This function is designed for summary statistics containing estimated regression coefficients and their estimated covariance matrices.
#'
#'
#' @details
#' Run a Gibbs sampler using a hierarchical spike.
#'
#'
#' @param Betah1 A numerical vector of length mg representing the regression coefficients for the first trait.
#' @param Betah2 A numerical vector of length mg representing the regression coefficients for the second trait.
#' @param Sigmah1 A mg*mg positive definite covariance matrix where it is estimated covariance matrix Betah1.
#' @param Sigmah2 A mg*mg Fpositive definite covariance matrix where it is estimated covariance matrix Betah2.
#' @param kappa0 Initial value for kappa (its dimension is equal to nchains).
#' @param kappastar0 Initial value for kappastar (its dimension is equal to nchains).
#' @param sigma20 Initial value for sigma2 (its dimension is equal to nchains).
#' @param s20 Initial value for s2 (its dimension is equal to nchains).
#' @param mg Number of variables in the group.
#' @param niter Number of iterations for the Gibbs sampler.
#' @param burnin Number of burn-in iterations.
#' @param nthin The lag of the iterations used for the posterior analysis is defined (or thinning rate).
#' @param nchains Number of Markov chains, when nchains>1, the function calculates the Gelman-Rubin convergence statistic, as modified by Brooks and Gelman (1998).
#' @param a1,a2 Hyperparameters of kappa. Default is a1=0.1 and a2=0.1.
#' @param c1,c2 Hyperparameters of kappastar. Default is c1=1 and c2=1.
#' @param d1,d2 Hyperparameters of sigma2. Default is d1=0.1 and d2=0.1.
#' @param e2 Initial value for doing Monte Carlo EM algorithm to estimate hyperparameter of s2.
#' @param snpnames Names of variables for the group.
#' @param genename Name of group.
#'
#' @return
#' - mcmcchain: The list of simulation output for all parameters.
#' - Criteria: genename, snpnames, the number of studies with nonzero signal by median and the number of studies with nonzero signal by credible interval (CI).
#' - Statistics of Trait 1 for Beta_1: Summary statistics including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#' - Statistics of Trait 2  for Beta_2: Summary statistics including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#' - Other Parameters: Summary statistics for kappa, kappastar and sigma2 including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#'
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' T. Baghfalaki, P.E. Sugier, T. Truong, A.T. Pettitt, K. Mengersen and B. Liquet. (2020). Bayesian meta-analysis models to Cross Cancer Genomic Investigation of pleiotropic effects using group structure. *Submitted in Statistics in Medicine*.
#'
#' @example inst/exampleHS.R
#'
#' @md
#' @export


HS=function(Betah1=Betah1, Betah2=Betah2, Sigmah1=Sigmah1, Sigmah2=Sigmah2, kappa0=kappa0, kappastar0=kappastar0, sigma20=sigma20, s20=s20,
            mg, niter=1000, burnin=500, nthin=2, nchains=2, a1=0.1, a2=0.1, d1=0.1, d2=0.1, c1=1, c2=1, e2=1, snpnames, genename){
  a1=a1; a2=a2; d1=d1; d2=d2; c1=c1; c2=c2; e2=e2
  Efinal=e2_Monte_Carlo_EM(Betah1=Betah1, Betah2=Betah2, Sigmah1=Sigmah1, Sigmah2=Sigmah2, kappa0 = kappa0[1],
                           kappastar0=kappastar0[1], sigma20 =sigma20[1],  s20 = s20[1],
                           mg=mg, a1=a1, a2=a2, d1=d1, d2=d2, c1=c1, c2=c2, e2=e2, snpnames, genename)

  RES1=list()
  RESBeta1= RESBeta2= matrix(0,1,mg)
  RESOthers=matrix(0,1,3)
  for(j in 1:nchains){

    RES1[[j]]=HS0(Betah1=Betah1, Betah2=Betah2, Sigmah1=Sigmah1, Sigmah2=Sigmah2, kappa0 = kappa0[j],
                  kappastar0=kappastar0[j], sigma20 = sigma20[j], s20=s20[j],
                  mg=mg,niter=niter, burnin=burnin, nthin=nthin, a1=a1, a2=a2, d1=d1, d2=d2, c1=c1, c2=c2,
                  e2=Efinal, snpnames=snpnames, genename=genename)


    RESBeta1=rbind(RESBeta1,RES1[[j]]$mcmcchain$Beta1)
    RESBeta2=rbind(RESBeta2,RES1[[j]]$mcmcchain$Beta2)
    RESOthers=rbind(RESOthers,cbind(RES1[[j]]$mcmcchain$kappa,RES1[[j]]$mcmcchain$kappastar,RES1[[j]]$mcmcchain$sigma2))

  }
  RESBeta1=RESBeta1[-1,]
  RESBeta2=RESBeta2[-1,]
  RESOthers=RESOthers[-1,]

  TabB1=wiqid::simpleRhat(RESBeta1,n.chains=nchains)
  TabB2=wiqid::simpleRhat(RESBeta2,n.chains=nchains)
  TabOthers=wiqid::simpleRhat(RESOthers,n.chains=nchains)

  TabOthers=t(TabOthers)
  colnames(TabOthers)=c("kappa","kappastar","sigma2")

  Tab=cbind(snpnames,TabB1,TabB2)
  colnames(Tab)<-c("Name of SNP", "BGR for Beta_1", "BGR for Beta_2")
  Tabb=list(Tab,TabOthers)


  RES1new=list(RES1,Tabb)
  names(RES1new) <- c("Outputs", "BGR" )

  ifelse(nchains==1,return(RES1),return(RES1new))

}



#' Internal: e2_Monte_Carlo_EM
#'
#' @param Betah1 A numerical vector of length mg representing the regression coefficients for the first trait.
#' @param Betah2 A numerical vector of length mg representing the regression coefficients for the second trait.
#' @param Sigmah1 A mg*mg positive definite covariance matrix where it is estimated covariance matrix Betah1.
#' @param Sigmah2 A mg*mg positive definite covariance matrix where it is estimated covariance matrix Betah2.
#' @param kappa0 Initial value for kappa.
#' @param kappastar0 Initial value for kappastar.
#' @param sigma20 Initial value for sigma2.
#' @param s20 Initial value for s2.
#' @param mg Number of variables in the group.
#' @param a1,a2 Hyperparameters of kappa. Default is a1=0.1 and a2=0.1.
#' @param c1,c2 Hyperparameters of kappastar. Default is c1=1 and c2=1.
#' @param d1,d2 Hyperparameters of sigma2. Default is d1=0.1 and d2=0.1.
#' @param e2 Initial value for doing Monte Carlo EM algorithm to estimate hyperparameter of s2.
#' @param snpnames Names of variables for the group.
#' @param genename Name of group.
#' @export
e2_Monte_Carlo_EM=function(Betah1=Betah1, Betah2=Betah2, Sigmah1=Sigmah1, Sigmah2=Sigmah2, kappa0 = kappa0,
                           kappastar0=kappastar0, sigma20 = sigma20, s20=s20,
                           mg, a1=a1, a2=a2, d1=d1, d2=d2, c1=c1, c2=c2, e2=e2, snpnames, genename){
  N00=100
  Beta1=b1=matrix(1,N00,mg)
  Beta2=b2=matrix(1,N00,mg)

  b1[1,]=Betah1
  b2[1,]=Betah2

  kappa=rep(0,N00);kappa[1]=kappa0
  kappastar=rep(0,N00);kappastar[1]=kappastar0

  sigma2=rep(1,N00);sigma2[1]=sigma20
  s2=rep(1,N00);s2[1]=s20

  tau1=tau2=matrix(1,N00,mg)

  MeanBeta1=MeanBeta2=MedianBeta1=MedianBeta2=c()

  Geneplotci=Geneplotmed=c()

  PROB1j=PROB2j=c()
  PROB1=PROB2=c()
  K=2
  for(r in 2:N00){
    V1=diag(tau1[r-1,]) ;  V2=diag(tau2[r-1,])

    ##################### Beta_k ####
    Omega1=V1%*%MASS::ginv(Sigmah1)%*%V1+(1/sigma2[r-1])*diag(mg)
    Mean1=MASS::ginv(Omega1)%*%V1%*%MASS::ginv(Sigmah1)%*%Betah1
    kappat1=(1/(1+(kappa[r-1]/(1-kappa[r-1])*exp(.5*t(Mean1)%*%Omega1%*%Mean1-mg/2*log(sigma2[r-1])-.5*determinant(Omega1, logarithm = TRUE)$modulus))))

    Tab=rbinom(1,1,as.numeric(kappat1))
    HH1=MASS::ginv(Omega1)
    gdata::lowerTriangle(HH1) <- gdata::upperTriangle(HH1, byrow=TRUE)
    b1[r,]=mvtnorm::rmvnorm(1, mean = Mean1, sigma =HH1)*(1-Tab)


    Omega2=V2%*%MASS::ginv(Sigmah2)%*%V2+(1/sigma2[r-1])*diag(mg)
    Mean2=MASS::ginv(Omega2)%*%V2%*%MASS::ginv(Sigmah2)%*%Betah2
    kappat2=(1/(1+(kappa[r-1]/(1-kappa[r-1])* exp(.5*t(Mean2)%*%Omega2%*%Mean2-mg/2*log(sigma2[r-1])-.5*determinant(Omega2, logarithm = TRUE)$modulus))))

    Tab=rbinom(1,1,as.numeric(kappat2))
    HH1=MASS::ginv(Omega2)
    gdata::lowerTriangle(HH1) <- gdata::upperTriangle(HH1, byrow=TRUE)
    b2[r,]=mvtnorm::rmvnorm(1, mean = Mean2, sigma =HH1)*(1-Tab)


    Beta1[r,]=b1[r,]%*%V1
    Beta2[r,]=b2[r,]%*%V2

    ##################### tau1 and tau2 ####
    m1=v1=ZZ1= m2=v2=ZZ2=c()
    sigmabar1=sigmabar2=c()
    for(j in 1:mg){
      sigmabar1[j]=Sigmah1[j,j]-Sigmah1[j,-j]%*%solve(Sigmah1[-j,-j])%*%Sigmah1[-j,j]
      vkj1=b1[r,j]^2/sigmabar1[j]+1/s2[r-1]
      A1=Betah1[j]-Sigmah1[j,-j]%*%solve(Sigmah1[-j,-j])%*%(Betah1[-j]-Beta1[r,-j])
      B1=A1*b1[r,j]/(vkj1*sigmabar1[j])
      vkj1B12=A1^2*b1[r,j]^2/(vkj1*sigmabar1[j]^2)
      T1=(log(2)+pnorm(B1*sqrt(vkj1), log.p = TRUE))*(B1<0)+0*(B1>=0)
      prob1=(1+kappastar[r-1]/(1-kappastar[r-1])*1/sqrt(vkj1*s2[r-1])*exp(vkj1B12/2+T1))^-1
      m1[j]=B1;v1[j]=1/vkj1;ZZ1[j]=rbinom(1,1,prob1)


      sigmabar2[j]=Sigmah2[j,j]-Sigmah2[j,-j]%*%solve(Sigmah2[-j,-j])%*%Sigmah2[-j,j]
      vkj2=b2[r,j]^2/sigmabar2[j]+1/s2[r-1]
      A2=Betah2[j]-Sigmah2[j,-j]%*%solve(Sigmah2[-j,-j])%*%(Betah2[-j]-Beta2[r,-j])
      B2=A2*b2[r,j]/(vkj2*sigmabar2[j])
      vkj2B22=A2^2*b2[r,j]^2/(vkj2*sigmabar2[j]^2)
      T2=(log(2)+pnorm(B2*sqrt(vkj2), log.p = TRUE))*(B2<0)+0*(B2>=0)
      prob2=(1+kappastar[r-1]/(1-kappastar[r-1])*1/sqrt(vkj2*s2[r-1])*exp(vkj2B22/2+T2))^-1
      m2[j]=B2;v2[j]=1/vkj2;ZZ2[j]=rbinom(1,1,prob2)
    }
    tau1[r,]=(1-ZZ1)*truncnorm::rtruncnorm(mg, a=0, b=Inf, mean = m1, sd = sqrt(v1))
    tau2[r,]=(1-ZZ2)*truncnorm::rtruncnorm(mg, a=0, b=Inf, mean = m2, sd = sqrt(v2))

    ##################### kappa ####
    K0=0
    if(max(b1[r,])==0)(K0=K0+1)
    if(max(b2[r,])==0)(K0=K0+1)
    kappa[r]=rbeta(1,K-K0+a1,K0+a2)
    #################### kappas ####
    Tai=c(tau1[r,],tau2[r,])
    index=1:length(Tai)
    L0=length(index[Tai==0])
    kappastar[r]=rbeta(1,K*mg-L0+c1,L0+c2)
    ##################### sigma2 ####
    sigma2[r]=invgamma::rinvgamma(1,shape=mg*(K-K0)/2+d1,rate=(Beta1[r,]%*%Beta1[r,]+Beta2[r,]%*%Beta2[r,])/2+d2)
    while(sigma2[r]>10){
      sigma2[r]=invgamma::rinvgamma(1,shape=mg*(K-K0)/2+d1,rate=(Beta1[r,]%*%Beta1[r,]+Beta2[r,]%*%Beta2[r,])/2+d2)
    }
    ##################### s2  ####
    s2[r]=invgamma::rinvgamma(1,shape=(K*mg-L0)/2+1,rate=(tau1[r,]%*%tau1[r,]+tau2[r,]%*%tau2[r,])/2+e2)
    e2=1/mean(1/s2[1:r])
  }
  ;e2
}



HS0=function(Betah1, Betah2, Sigmah1, Sigmah2, kappa0, kappastar0, sigma20, s20=s20, mg, niter=1000, burnin=500, nthin=2, a1=a1, a2=a2, c1=c1, c2=c2, d1=d1, d2=d2, e2=e2, snpnames, genename){

  Beta1=b1=matrix(0,niter,mg)
  Beta2=b2=matrix(0,niter,mg)

  b1[1,]=Betah1
  b2[1,]=Betah2

  kappa=rep(0,niter);kappa[1]=kappa0
  kappastar=rep(0,niter);kappastar[1]=kappastar0

  sigma2=rep(1,niter);sigma2[1]=sigma20
  s2=rep(1,niter);s2[1]=s20

  tau1=tau2=matrix(1,niter,mg)
  MeanBeta1=MeanBeta2=MedianBeta1=MedianBeta2=c()

  Geneplotci=Geneplotmed=c()

  PROB1j=PROB2j=c()
  PROB1=PROB2=c()
  K=2
  for(r in 2:niter){
    V1=diag(tau1[r-1,]) ;  V2=diag(tau2[r-1,])

    ##################### Beta_k ####
    Omega1=V1%*%MASS::ginv(Sigmah1)%*%V1+(1/sigma2[r-1])*diag(mg)
    Mean1=MASS::ginv(Omega1)%*%V1%*%MASS::ginv(Sigmah1)%*%Betah1
    kappat1=(1/(1+(kappa[r-1]/(1-kappa[r-1])*exp(.5*t(Mean1)%*%Omega1%*%Mean1-mg/2*log(sigma2[r-1])-.5*determinant(Omega1, logarithm = TRUE)$modulus))))
    Tab=rbinom(1,1,as.numeric(kappat1))
    HH1=MASS::ginv(Omega1)
    gdata::lowerTriangle(HH1) <- gdata::upperTriangle(HH1, byrow=TRUE)
    b1[r,]=mvtnorm::rmvnorm(1, mean = Mean1, sigma =HH1)*(1-Tab)


    Omega2=V2%*%MASS::ginv(Sigmah2)%*%V2+(1/sigma2[r-1])*diag(mg)
    Mean2=MASS::ginv(Omega2)%*%V2%*%MASS::ginv(Sigmah2)%*%Betah2
    kappat2=(1/(1+(kappa[r-1]/(1-kappa[r-1])* exp(.5*t(Mean2)%*%Omega2%*%Mean2-mg/2*log(sigma2[r-1])-.5*determinant(Omega2, logarithm = TRUE)$modulus))))
    Tab=rbinom(1,1,as.numeric(kappat2))
    HH1=MASS::ginv(Omega2)
    gdata::lowerTriangle(HH1) <- gdata::upperTriangle(HH1, byrow=TRUE)
    b2[r,]=mvtnorm::rmvnorm(1, mean = Mean2, sigma =HH1)*(1-Tab)


    Beta1[r,]=b1[r,]%*%V1
    Beta2[r,]=b2[r,]%*%V2

    ##################### tau1 and tau2 ####
    m1=v1=ZZ1= m2=v2=ZZ2=c()
    sigmabar1=sigmabar2=c()
    for(j in 1:mg){
      sigmabar1[j]=Sigmah1[j,j]-Sigmah1[j,-j]%*%solve(Sigmah1[-j,-j])%*%Sigmah1[-j,j]
      vkj1=b1[r,j]^2/sigmabar1[j]+1/s2[r-1]
      A1=Betah1[j]-Sigmah1[j,-j]%*%solve(Sigmah1[-j,-j])%*%(Betah1[-j]-Beta1[r,-j])
      B1=A1*b1[r,j]/(vkj1*sigmabar1[j])
      vkj1B12=A1^2*b1[r,j]^2/(vkj1*sigmabar1[j]^2)
      T1=(log(2)+pnorm(B1*sqrt(vkj1), log.p = TRUE))*(B1<0)+0*(B1>=0)
      prob1=(1+kappastar[r-1]/(1-kappastar[r-1])*1/sqrt(vkj1*s2[r-1])*exp(vkj1B12/2+T1))^-1
      m1[j]=B1;v1[j]=1/vkj1;ZZ1[j]=rbinom(1,1,prob1)


      sigmabar2[j]=Sigmah2[j,j]-Sigmah2[j,-j]%*%solve(Sigmah2[-j,-j])%*%Sigmah2[-j,j]
      vkj2=b2[r,j]^2/sigmabar2[j]+1/s2[r-1]
      A2=Betah2[j]-Sigmah2[j,-j]%*%solve(Sigmah2[-j,-j])%*%(Betah2[-j]-Beta2[r,-j])
      B2=A2*b2[r,j]/(vkj2*sigmabar2[j])
      vkj2B22=A2^2*b2[r,j]^2/(vkj2*sigmabar2[j]^2)
      T2=(log(2)+pnorm(B2*sqrt(vkj2), log.p = TRUE))
      prob2=(1+kappastar[r-1]/(1-kappastar[r-1])*1/sqrt(vkj2*s2[r-1])*exp(vkj2B22/2+T2))^-1
      m2[j]=B2;v2[j]=1/vkj2;ZZ2[j]=rbinom(1,1,prob2)
    }
    tau1[r,]=(1-ZZ1)*truncnorm::rtruncnorm(mg, a=0, b=Inf, mean = m1, sd = sqrt(v1))
    tau2[r,]=(1-ZZ2)*truncnorm::rtruncnorm(mg, a=0, b=Inf, mean = m2, sd = sqrt(v2))


    ##################### kappa ####
    K0=0
    if(max(b1[r,])==0)(K0=K0+1)
    if(max(b2[r,])==0)(K0=K0+1)

    kappa[r]=rbeta(1,K-K0+a1,K0+a2)
    ##################### kappas ####
    Tai=c(tau1[r,],tau2[r,])
    index=1:length(Tai)
    L0=length(index[Tai==0])
    kappastar[r]=rbeta(1,K*mg-L0+c1,L0+c2)
    ##################### sigma2 ####
    sigma2[r]=invgamma::rinvgamma(1,shape=mg*(K-K0)/2+d1,rate=(Beta1[r,]%*%Beta1[r,]+Beta2[r,]%*%Beta2[r,])/2+d2)
    while(sigma2[r]>10){
      sigma2[r]=invgamma::rinvgamma(1,shape=mg*(K-K0)/2+d1,rate=(Beta1[r,]%*%Beta1[r,]+Beta2[r,]%*%Beta2[r,])/2+d2)
    }
    ##################### s2  ####
    s2[r]=invgamma::rinvgamma(1,shape=(K*mg-L0)/2+1,rate=(tau1[r,]%*%tau1[r,]+tau2[r,]%*%tau2[r,])/2+e2)
  }
  #$$$$$$$$$$$$$$$$$$$$$$$$$ New Postorior samples After Burn-in $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  niter1=burnin+1
  index=niter1:niter
  indexn=index[(index %% nthin) == 0]

  kappastar=kappastar[indexn]
  kappa=kappa[indexn]
  sigma2=sigma2[indexn]
  Beta1=Beta1[indexn,]
  Beta2=Beta2[indexn,]
  Tau1=tau1[indexn,]
  Tau2=tau2[indexn,]
  s2=s2[indexn]


  MeanBeta1=apply(Beta1,2,mean)
  MeanBeta2=apply(Beta2,2,mean)

  SDBeta1=apply(Beta1,2,sd)
  SDBeta2=apply(Beta2,2,sd)

  QBeta1=t(apply(Beta1, 2, function(x) quantile(x, c(.025, 0.5, .975))))
  QBeta2=t(apply(Beta2, 2, function(x) quantile(x, c(.025, 0.5, .975))))


  Other=cbind(kappa,kappastar,sigma2,s2)
  MeanOther=apply(Other,2,mean)

  SDOther=apply(Other,2,sd)

  QOther=t(apply(Other, 2, function(x) quantile(x, c(.025, 0.5, .975))))



  PGENE1m=PGENE2m=rep(0,mg)
  B1CI=B2CI=c()
  for(int in 1:mg){
    B1CI=median(Beta1[,int])
    B2CI=median(Beta2[,int])

    if(B1CI!=0)(PGENE1m[int]=1)
    if(B2CI!=0)(PGENE2m[int]=1)
  }

  Geneplotmed=PGENE1m+PGENE2m
  ###########################################
  PGENE1=PGENE2=rep(0,mg)
  B1CI=B2CI=c()
  for(int in 1:mg){
    B1CI=quantile(Beta1[,int],c(0.025,0.975))
    B2CI=quantile(Beta2[,int],c(0.025,0.975))
    if((0< B1CI[1]) | (0>B1CI[2]))(PGENE1[int]=1)
    if((0<B2CI[1]) | (0>B2CI[2]))(PGENE2[int]=1)
  }

  Geneplotci=PGENE1+PGENE2

  mcmcchain=list(kappa=kappa,kappastar=kappastar,sigma2=sigma2,Beta1=Beta1,Beta2=Beta2,tau1=Tau1,tau2=Tau2,s2=s2)


  Reslast= list("Name of Gene"=genename,"Name of SNP"=snpnames,
                "# studies nonzero signal by CI"=Geneplotci,"# studies nonzero signal by Med"=Geneplotmed)


  Trait_1= cbind(snpnames, MeanBeta1, SDBeta1, QBeta1)
  colnames(Trait_1)= cbind("Name of SNP", "Mean", "SD", "val2.5pc",	"Median",	"val97.5pc")

  Trait_2= cbind(snpnames, MeanBeta2, SDBeta2, QBeta2)
  colnames(Trait_2)= cbind("Name of SNP", "Mean", "SD", "val2.5pc",	"Median",	"val97.5pc")


  Others= cbind(MeanOther, SDOther, QOther)
  colnames(Others)= cbind( "Mean", "SD", "val2.5pc",	"Median",	"val97.5pc")

  OUT=list(sim.matrix=mcmcchain, Criteria=Reslast, Trait_1, Trait_2, Others)

  names(OUT) <- c("mcmcchain", "Criteria", "Statistics of Trait 1 for Beta_1", "Statistics of Trait 2  for Beta_2", "Other Parameters")

  ;OUT
}


