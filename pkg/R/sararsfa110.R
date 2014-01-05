sararsfa110.params = function(parameters){
  k = length(parameters)
  beta <- parameters[1:(k-4)]
  sigmaV <- parameters[k-3]
  sigmaU <- parameters[k-2]
  rho <- parameters[k-1]
  rho2 <- parameters[k]
  
  names(beta) = paste("Beta", seq(k-4), sep = "")
  names(sigmaV) = "SigmaV"
  names(sigmaU) = "SigmaU"
  names(rho) = "Rho"
  names(rho2) = "Rho2"
  
  return(list(beta=beta,sigmaV=sigmaV,sigmaU=sigmaU,rho=rho,rho2=rho2))
}

sararsfa110.hnormal.lf <- function(parameters){ 
  y = envir.get("y")
  X = envir.get("X")
  W = envir.get("W")
  W2 = envir.get("W2")
  p = sararsfa110.params(parameters)
  counter = envir.counter("sararsfa110.hnormal.lf")
  logger.debug(paste("Evaluating 'sararsfa110.hnormal.lf' for parameters","[run=",counter,"]:"),parameters) 
  ret = -10e8
  if ((p$sigmaV>0) && (p$sigmaU>0) && (abs(p$rho)<1) && (abs(p$rho2)<1)){ 
    n = length(y)
    e <- y  - p$rho * W %*% y -  X %*% p$beta
    I <- diag(n)    
    Sp2 <- solve(I-p$rho2*W2)
    mCapitalSigma <- t(Sp2)%*%Sp2
    mCapitalTheta <- p$sigmaU^2*I+p$sigmaV^2*mCapitalSigma
    mCapitalOmega <-p$sigmaU^2*p$sigmaV^2*mCapitalSigma %*% solve(mCapitalTheta)
    mMu <- -p$sigmaV^(-2)*mCapitalOmega%*%solve(mCapitalSigma)%*%e
    ret<- log(det(I-p$rho*W))+log(ptmvnorm(lowerx=rep(0, n),upperx = rep( Inf, n),mean=as.vector(t(mMu)), sigma=mCapitalOmega))+log(dmvnorm(x=as.vector(e),mean=rep(0, n), sigma=mCapitalTheta))
  }else{
   logger.debug("Parameters are out of space")
  }
  if (is.nan(ret) || (ret==-Inf)) ret = -10e8
  logger.debug(paste("'sararsfa110.hnormal.lf' value:",-ret))
  return(-ret)
}


sararsfa110.hnormal.ini<-function(formula, data){
  logger.debug("Calculating initial values")
  W = envir.get("W")
  sarsfa <- sarsfa.hnormal.estimator(formula, data, W)
  
  if (length(sarsfa$beta)==0) {
    print("sarsfa is failed")
    envir.assign("est.failed", TRUE)
    return(NULL)
  }
  rho2 = 0
  result <-c(sarsfa$beta,sarsfa$sigmaV,sarsfa$sigmaU,sarsfa$rho, rho2)
  names(result) <-c(paste("Beta", seq(from=0, to=length(sarsfa$beta)-1),sep = ""),"SigmaV","SigmaU", "Rho", "Rho2")
  logger.info("Initial values:", result)
  return(result)
}

sararsfa110.hnormal.estimator <- function(formula, data,W,W2,logging = "quiet", ini.values = NULL,control=NULL){
  envir.init()
  envir.assign("logging.level",logging)
  if (!is.null(control)) envir.assign("control",control)
  logger.start()
  logger.debug("Estimator started")
  
  sfa.prepare(formula, data)
  envir.assign("W", W)
  envir.assign("W2", W2)
  
  if (is.null(ini.values)){
    ini.values = sararsfa110.hnormal.ini(formula, data)  
  } 
  
  estimates <- optim.estimator(formula, data, sararsfa110.hnormal.lf, ini.values, gr=NULL)
  
  est_failed = !is.null(envir.get("est.failed"))
  res = sararsfa110.params(estimates$estimate)
  if(est_failed){
    res$failed = T
  }
  logger.info("Estimates:",c(res$beta,res$sigmaV,res$sigmaU,res$rho,res$rho2))  
  print("")
  envir.finalise()
  return(res)
}