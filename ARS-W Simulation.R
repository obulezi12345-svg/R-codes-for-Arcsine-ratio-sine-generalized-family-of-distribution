## Clean up everything
remove(list=objects())

## package
library(BMT)
library(bbmle)

## pdf of NoCS-W distribution
fx = function(x, a, t, l){(((pi)*a*cos((pi/2)*(1-exp(-t*x^(l)))))*t*l*x^(l-1)*exp(-t*x^l)/(a+(2-a)*sin((pi/2)*(1-exp(-t*x^(l)))))^2)}
## cdf of the NOCS-W distribution
Fx = function(x, a, t, l){1-((a*(1-sin((pi/2)*(1-exp(-t*x^(l))))))/(a+(2-a)*sin((pi/2)*(1-exp(-t*x^(l))))))}

## survival function of OLiP distribution
Sx = function(x, a, t, l){1-Fx(x,a, t, l)}

## quantile (Non-closed form)
## F(x)-u = 0
xp = function(u,a, t, l){
  ## cdf with x
  F1 = function(x){Fx(x,a, t, l)}
  ## inverse cdf 
  Inv = function(u){uniroot(function(x){F1(x) - u}, c(0,1),
                            extendInt = "yes")$root}
  ## generating data 
  h = c()
  for(j in 1:length(u)){h[j] = Inv(u[j])}
  return(h)
}

## =======================================================================
##                       Method of estimation functions                 ##
## =======================================================================
## 1- Negative log-likehood function (MLE) -------------------------------
Nlog_like = function(parm){
  a = parm[1]
  t = parm[2]
  l = parm[3]
  Nlog_like = -sum(log(fx(x,a, t, l)))
  return(Nlog_like)
}

## 2- Maximum product of spacing (MPS) estimation ------------
dNOCSW <- function(x,a, t, l){fx(x,a, t, l)}
pNOCSW <- function(q,a, t, l){Fx(q,a, t, l)}
qNOCSW <- function(p,a, t, l){qx(p,a, t, l)}

## 3- Least square (LS) -------------------------------------
LS = function(parm){
  a = parm[1]
  t = parm[2]
  l  = parm[3]
  LS = sum((Fx(x,a, t, l)-((1:n)/(n+1)))^2)
  return(LS)
}

## 4- Weighted least square method (WLS) ----------------------
WLS = function(parm){
  a = parm[1]
  t = parm[2]
  l  = parm[3]
  j = 1:n
  w = (n+1)^2*(n+2)/(j*(n-j+1))
  WLS = sum(w*(Fx(x,a, t, l)-((j)/(n+1)))^2)
  return(WLS)
}

## 5- Method of Cramer - Von Mises (CVM) -----------------------
CVM = function(parm){
  a = parm[1]
  t = parm[2]
  l  = parm[3]
  j = 1:n
  w = 1/(12*n)
  CVM = w + sum((Fx(x,a, t, l)-((2*j-1)/(2*n)))^2)
  return(CVM)
}

## 6- Anderson-Darling estimation (AD) --------------------------
AD = function(parm){
  a = parm[1]
  t     = parm[2]
  l  = parm[3]
  j = 1:length(x)
  A = log(Fx(x[j],a, t, l))+log(Sx(x[n+1-j],a, t, l))
  AD = -n-(1/n)*sum((2*j-1)*A)
  return(AD)
}


## 7- Right Tail Anderson-Darling estimation (RTAD) ---------------
RTAD = function(parm){
  a = parm[1]
  t = parm[2]
  l  = parm[3]
  j = 1:length(x)
  R1 = 2*sum(Fx(x[j],a, t, l)) 
  R2 = log(Sx(x[n+1-j],a, t, l))
  RTAD = (n/2)-R1-(1/n)*sum((2*j-1)*R2)
  return(RTAD)
}

## 8- Percentile estimators (PS) ----------------------------------
PERC = function(parm){
  a = parm[1]
  t = parm[2]
  l  = parm[3]
  j = 1:n
  u = j/(n+1)
  PERC = sum((x-xp2(u,a, t, l))^2, na.rm = TRUE)
  return(PERC)
}


## function of resluts
Rslu = function(Estimate,initial){
  Mean = mean(Estimate, na.rm = TRUE)
  MSE  = mean((Estimate-initial)^2, na.rm = TRUE)
  RMSE = sqrt(MSE)
  Bais = abs(Mean-initial)
  Rslu = c(initial, Mean, MSE, RMSE,Bais)
}


## =======================================================================
##                          Parmaeters                                  ##
## =======================================================================
## sample size 
N = c(50, 100, 150 , 200)
#N=100
## parameters of SAWe dist.
Parm11 = c(20.6,1.06,2,6)
Parm12 = c(9.2,0.5,0.2)
Parm13 = c(4.1,0.1,6.2)
Parm14 = c(8.2,1.2,0,6)
Parm21 = c(0.5,0.02,1.27)
Parm22 = c(0.01,4,0.001)
Parm23 = c(0.7,0.35,4)
Parm24 = c(0.2,1.09,5.02)
#Parmnew = c(0.50,0.75,1.50)

Parm = Parm11
a  = Parm[1]
t = Parm[2]
l = Parm[3]

Case = "Parm11"


## =======================================================================
##                               Simulation                             ##
## =======================================================================
Mat = matrix(NA, nrow =  24, ncol = 12)
Rl = 1

for(kk in 1:length(N)){
  n = N[kk]
  
  ## N. of simualtion
  N.sim = 1000
  
  ## saving vectors
  MLE   = matrix(NA, nrow = N.sim, ncol = 3)
  LSE   = matrix(NA, nrow = N.sim, ncol = 3)
  WLSE  = matrix(NA, nrow = N.sim, ncol = 3)
  MPSE  = matrix(NA, nrow = N.sim, ncol = 3)
  CVME  = matrix(NA, nrow = N.sim, ncol = 3)
  ADE   = matrix(NA, nrow = N.sim, ncol = 3)
  RTADE = matrix(NA, nrow = N.sim, ncol = 3)
  PERCE = matrix(NA, nrow = N.sim, ncol = 3)
  
  
  ## loopping of simualtion
  for(i in 1:N.sim){
    ## for fixed random generating 
    #set.seed(i)
    
    ## generating uinform
    u = sort(runif(n))
    
    ## generate data from OLiP dist
    x = xp(u,a, t, l)
    
    ## Estimate of theta for MLE and MPS
    #MLE[i,4] = min(x)
    
    ## Maximum likehood estimation 
    MLE[i,1:3] =  suppressWarnings(optim(par=c(a=a, t=t, l=l), fn = Nlog_like, 
                                         hessian = TRUE, method = "N")$par)
    
    ## Maximum Product Spacing estimation
    tryCatch({
      Fit_mps  = mpsedist(data =  x, distr = "NOCSW", start = list(a=a, t=t, l=l),optim.method = "N")
      MPSE[i,] = as.numeric(Fit_mps$estimate)},
      error = function(e) {cat("ERROR :",conditionMessage(e), "\n")})
    
    ## Least square estimation
    LSE[i,]  = suppressWarnings(optim(par=c(a=a, t=t, l=l), fn = LS, method = "N")$par)
    
    ## Weighted Least square estimation
    WLSE[i,] = suppressWarnings(optim(par=c(a=a, t=t, l=l), fn = WLS, method = "N")$par)
    
    ## Cramer von Mises estimation
    CVME[i,] = suppressWarnings(optim(par=c(a=a, t=t, l=l), fn = CVM, method = "N")$par)
    
    ## Anderson- Darling estimation
    ADE[i,] = suppressWarnings(optim(par=c(a=a, t=t, l=l), fn = AD, method = "N")$par)
    
    ## Right Tail Anderson- Darling estimation
    RTADE[i,] = suppressWarnings(optim(par=c(a=a, t=t, l=l), fn = RTAD, method = "N")$par)
    
    # ## Percentile Estimators
    # PERCE[i,] = suppressWarnings(optim(par=c(lambda=lambda,k=k,theta=theta), fn = PERC, method = "N")$par)
  }
  
  
  ##==================================================================================
  ##                                Print Resluts                                   ##
  ##==================================================================================
  #1# estimate uinsg Monte-Carlo: MLE
  a_mle = Rslu(MLE[,1],a)
  t_mle      = Rslu(MLE[,2],t)
  l_mle  = Rslu(MLE[,3],l)
  Name = c("initial","Mean", "MSE","RMSE" ,"Bias")
  ReslutofMLE = data.frame(Name,a_mle,t_mle,l_mle)
  
  #2# estimate uinsg Monte-Carlo: MPSE
  a_mpse = Rslu(MPSE[,1],a)
  t_mpse      = Rslu(MPSE[,2],t)
  l_mpse  = Rslu(MPSE[,3],l)
  ReslutofMPS = data.frame(Name,a_mpse,t_mpse,l_mpse)
  
  #3# estimate uinsg Monte-Carlo: LSE
  a_ls = Rslu(LSE[,1],a)
  t_ls      = Rslu(LSE[,2],t)
  l_ls  = Rslu(LSE[,3],l)
  ReslutofLS = data.frame(Name,a_ls,t_ls,l_ls)
  
  #4# estimate uinsg Monte-Carlo: WLSE
  a_wlse  = Rslu(WLSE[,1],a)
  t_wlse       = Rslu(WLSE[,2],t)
  l_wlse   = Rslu(WLSE[,3],l)
  ReslutofWLSE = data.frame(Name,a_wlse,t_wlse,l_wlse)
  
  #5# estimate uinsg Monte-Carlo: CVME
  a_cvme = Rslu(CVME[,1],a)
  t_cvme      = Rslu(CVME[,2],t)
  l_cvme  = Rslu(CVME[,3],l)
  ReslutofCVM = data.frame(Name,a_cvme,t_cvme,l_cvme)
  
  #6# estimate uinsg Monte-Carlo: ADE
  a_ade = Rslu(ADE[,1],a)
  t_ade      = Rslu(ADE[,2],t)
  l_ade  = Rslu(ADE[,3],l)
  ReslutofADE = data.frame(Name,a_ade,t_ade,l_ade)
  
  #7# estimate uinsg Monte-Carlo: RTADE
  a_rtade = Rslu(RTADE[,1],a)
  t_rtade      = Rslu(RTADE[,2],t)
  l_rtade  = Rslu(RTADE[,3],l)
  ReslutofRTADE = data.frame(Name,a_rtade,t_rtade,l_rtade)
  
  # #8# estimate uinsg Monte-Carlo: PERCE
  # lambda_ps = Rslu(PERCE[,1],lambda)
  # k_ps      = Rslu(PERCE[,2],k)
  # theta_ps  = Rslu(PERCE[,3],theta)
  # ReslutofPE  = data.frame(Name,lambda_ps,k_ps,theta_ps)
  
  
  ## Summary of Estimation
  Descrption = c("Sample size","N.simluation","Initial a", "Initial t","Initial l")
  value      = c(n,N.sim,a,t,l)
  Summary    = data.frame(Descrption,value)
  
  ## Print resluts
  print(Summary)
  print(ReslutofMLE)
  print(ReslutofMPS)
  print(ReslutofLS)
  print(ReslutofWLSE)
  print(ReslutofCVM)
  print(ReslutofADE)
  print(ReslutofRTADE)
  
  
  for(gg in 2:4){
    Mat[(gg-1), Rl:(Rl+1)] = c(ReslutofMLE[c(5,5),gg])
    Mat[(gg+3), Rl:(Rl+1)] = c(ReslutofMPS[c(5,5),gg])
    Mat[(gg+7), Rl:(Rl+1)] = c(ReslutofLS[c(5,5),gg])
    Mat[(gg+11),Rl:(Rl+1)] = c(ReslutofWLSE[c(5,5),gg])
    Mat[(gg+15),Rl:(Rl+1)] = c(ReslutofCVM[c(5,5),gg])
    Mat[(gg+19),Rl:(Rl+1)] = c(ReslutofADE[c(5,5),gg])
    Mat[(gg+23),Rl:(Rl+1)] = c(ReslutofRTADE[c(5,5),gg])
  }
  
  Rl = Rl+3
  print(Mat)
}


write.csv(Mat,file=Case)
####### in tables we are put the Mean square error and 

