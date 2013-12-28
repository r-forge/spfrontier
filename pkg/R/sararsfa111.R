sararsfa111.tnormal.lf <- function(parameters,env){
  num = 0
  if (exists("sararsfa111.tnormal.lf_number", envir =env)){
    num = get("sararsfa111.tnormal.lf_number", envir =env)
  }
  num = num +1 
  assign("sararsfa111.tnormal.lf_number", num, envir=env)
  print(paste("Lf number ",num))
  
  y = get("y", envir =env)
  X = get("X", envir =env)
  W = get("W", envir =env)
  W2 = get("W2", envir =env)
  W3 = get("W3", envir =env)
  n = length(y)
  maxValue= get("max_ll_value", envir =env)
  k = length(parameters)
  
  beta <- parameters[1:(k-6)]
  sigmaV <- parameters[k-5]
  sigmaU <- parameters[k-4]
  
  
  rho <- parameters[k-3]
  rho2 <- parameters[k-2]
  rho3 <- parameters[k-1]
  mu <- parameters[k]
  
  
  if (sigmaV<=0) return (maxValue)
  if (sigmaU<=0) return (maxValue)
  if (abs(rho)>1) return (maxValue)
  if (abs(rho2)>1) return (maxValue)
  if (abs(rho3)>1) return (maxValue)
  
  e <- y  - rho * W %*% y -  X %*% beta
  I <- diag(n)
  
  Sp2 <- solve(I-rho2*W2)
  Sp3 <- solve(I-rho3*W3)
  
  mSigma = sigmaV^2*t(Sp2)%*%Sp2
  mOmega = sigmaU^2*t(Sp3)%*%Sp3
  mC = mSigma + mOmega
  mB = mOmega%*%solve(mC)%*%mSigma
  #mB = mOmega%*%chol2inv(chol(mC))%*%mSigma
  
  mA = mB %*% solve(mSigma)
  mD = -mB %*% solve(mOmega)
  vMu = rep(mu, n)
  
  logl <- maxValue
  tryCatch({
  logl<- - (log(det(I-rho*W))-log(ptmvnorm(lowerx=rep(-Inf, n),upperx = rep(0, n),mean=as.vector(t(vMu)), sigma=mOmega))
            +log(dmvnorm(x=as.vector(e+vMu),mean=rep(0, n), sigma=mC))
            +log(ptmvnorm(lowerx=rep(-Inf, n),upperx = rep(0, n),mean=as.vector(t(-mA%*%e-mD%*%vMu)), sigma=mB))
            )
  }, error = function(e){
    print(e$message)
  })
  if (is.nan(logl) || logl>maxValue) {
    logl <-maxValue
  }
  names(logl) <- "Log-Lik sararsfa111 normal/truncated-normal"
  return(logl)
}


sararsfa111.tnormal.ini<-function(formula, data, env=new.env()){
  
  W = get("W", envir =env)
  
  sarsfa <- sarsfa.hnormal.estimator(formula, data, W, env=env)
  
  
  sigma <- sarsfa$sigma
  rm("est_success", envir=env)
  if (length(sigma)==0) {
    print("sarsfa is failed")
    assign("est_success", FALSE, envir=env)
    return(NULL)
  }
  
  rho2 = 0.2
  rho3 = 0.3
  mu = 1
  result <-c(sarsfa$beta,sarsfa$sigmaV,sarsfa$sigmaU,sarsfa$rho, rho2, rho3, mu)
  names(result) <-c(paste("Beta", seq(from=0, to=length(gamma)+1),sep = ""),"SigmaV","SigmaU", "Rho", "Rho2", "Rho3", "Mu")
  
  return(result)
}

sararsfa111.tnormal.estimator <- function(formula, data,W,W2,W3,silent=TRUE,env=new.env()){
  assign("max_ll_value", 1E10, envir=env)
  assign("W", W, envir=env)
  assign("W2", W2, envir=env)
  assign("W3", W3, envir=env)
  num = 0
  if (exists("seq_number", envir =env)){
    num = get("seq_number", envir =env)
  }
  num = num +1 
  assign("seq_number", num, envir=env)
  print(paste("Run number ",num))
  estimates <- optim.estimator(formula, data, sararsfa111.tnormal.lf, sararsfa111.tnormal.ini, gr=NULL, silent=silent,env=env)
  if(!is.null(estimates)){
    est <-as.vector(estimates$estimate)
    k <- length(est)
    beta <- as.vector(est[1:(k-6)])
    sigmaV <- as.numeric(est[k-5])
    sigmaU <-  as.numeric(est[k-4])
    rho <-  as.numeric(est[k-3])
    rho2 <-  as.numeric(est[k-2])
    rho3 <-  as.numeric(est[k-1])
    mu <-  as.numeric(est[k])
    res <- list(beta=beta,rho=rho,rho2=rho2,rho3=rho3,mu=mu,sigmaV=sigmaV, sigmaU=sigmaU)
    return(res)
  } else{
    return(list())
  }
}