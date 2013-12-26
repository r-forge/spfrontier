sararsfa110.hnormal.lf <- function(parameters,env){
  y = get("y", envir =env)
  X = get("X", envir =env)
  W = get("W", envir =env)
  W2 = get("W2", envir =env)
  
  maxValue= get("max_ll_value", envir =env)
  k = length(parameters)
  
  beta <- parameters[1:(k-4)]
  sigmaV <- parameters[k-3]
  sigmaU <- parameters[k-2]
  
  
  rho <- parameters[k-1]
  rho2 <- parameters[k]
  
  e <- y  - rho * W %*% y -  X %*% beta
  I <- diag(n)
  
  Sp2 <- solve(I-rho2*W2)
  
  mCapitalSigma <- t(Sp2)%*%Sp2
  mCapitalTheta <- sigmaU^2*I+sigmaV^2*mCapitalSigma
  mCapitalOmega <-sigmaU^2*sigmaV^2*mCapitalSigma %*% solve(mCapitalTheta)
  mMu <- -sigmaV^(-2)*mCapitalOmega%*%solve(mCapitalSigma)%*%e
  logl<- - (log(det(I-rho*W))+log(pmvnorm(lower=rep(0, n),mean=as.vector(t(mMu)), sigma=mCapitalOmega))+log(dmvnorm(x=as.vector(e),mean=rep(0, n), sigma=mCapitalTheta)))
  
  if (is.nan(logl) || logl>maxValue) {
    logl <-maxValue
  }
  names(logl) <- "Log-Lik SFA normal/half-normal"
  return(logl)
}


sararsfa110.hnormal.ini<-function(formula, data, env=new.env()){
  
  W = get("W", envir =env)
  
  sarsfa <- sarsfa.hnormal.estimator(formula, data, W, env=env)
  
  
  sigma <- sarsfa$sigma
  rm("est_success", envir=env)
  if (length(sigma)==0) {
    print("sarsfa is failed")
    assign("est_success", FALSE, envir=env)
    return(NULL)
  }
  
  rho2 = 0
  result <-c(sarsfa$beta,sarsfa$sigmaV,sarsfa$sigmaU,sarsfa$rho, rho2)
  names(result) <-c(paste("Beta", seq(from=0, to=length(gamma)+1),sep = ""),"SigmaV","SigmaU", "Rho", "Rho2")
  
  return(result)
}

sararsfa110.hnormal.estimator <- function(formula, data,W,W2,silent=TRUE,env=new.env()){
  assign("max_ll_value", 1E10, envir=env)
  assign("W", W, envir=env)
  assign("W2", W2, envir=env)
  num = 0
  if (exists("seq_number", envir =env)){
    num = get("seq_number", envir =env)
  }
  num = num +1 
  assign("seq_number", num, envir=env)
  print(paste("Run number ",num))
  estimates <- optim.estimator(formula, data, sararsfa110.hnormal.lf, sararsfa110.hnormal.ini, gr=NULL, silent=silent,env=env)
  if(!is.null(estimates)){
    est <-as.vector(estimates$estimate)
  k <- length(est)
  beta <- as.vector(est[1:(k-4)])
  sigmaV <- as.numeric(est[k-3])
  sigmaU <-  as.numeric(est[k-2])
  rho <-  as.numeric(est[k-1])
  rho2 <-  as.numeric(est[k])
  res <- list(beta=beta,rho=rho,rho2=rho2,sigmaV=sigmaV, sigmaU=sigmaU)
  return(res)
  } else{
    return(list())
  }
}