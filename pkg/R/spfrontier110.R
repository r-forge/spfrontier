spfrontier110HN.params = function(parameters){
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

spfrontier110HN.logL <- function(parameters){ 
  counter = envir.counter("spfrontier110HN.logL")
  logger.debug(paste("Evaluating 'spfrontier110HN.logL' for parameters","[run=",counter,"]:"),parameters) 
  y = envir.get("y")
  X = envir.get("X")
  W = envir.get("W")
  W2 = envir.get("W2")
  p = spfrontier110HN.params(parameters)
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
  logger.debug(paste("spfrontier110HN.logL =",-ret))
  return(-ret)
}


spfrontier110HN.ini<-function(formula, data){
  logger.debug("spfrontier110HN: calculating initial values")
  W = envir.get("W")
  iniModel <- spfrontier(formula,data,model="spfrontier100HN",W_y=W)
  coefs <- iniModel@coefficients
  if (status(iniModel) > 0){
    iniModel <- spfrontier(formula,data,model="frontierHN")
    coefs <- iniModel@coefficients
    coefs$rho <- 0
  }
  if (length(coefs$beta)==0) {
    print("coefs is failed")
    envir.assign("est.failed", TRUE)
    return(NULL)
  }
  rho2 = 0
  result <-c(coefs$beta,coefs$sigmaV,coefs$sigmaU,coefs$rho, rho2)
  names(result) <-c(paste("Beta", seq(from=0, to=length(coefs$beta)-1),sep = ""),"SigmaV","SigmaU", "Rho", "Rho2")
  logger.info("Initial values:", result)
  return(result)
}

spfrontier110HN.prepare = function(formula, data,W,W2,...){
  prepareXY(formula, data)
  envir.assign("W", W)
  envir.assign("W2", W2)
}

spfrontier110HN.handle.estimates <- function(estimates){
  logger.debug("spfrontier110HN: Handling estimates")
  coefs <- list()
  logL <- 0
  if (!is.null(estimates$estimate)){
    status <- 0
    coefs <- spfrontier110HN.params(as.vector(estimates$estimate))
    
    logL <- estimates$value
    logger.info("Estimates:",unlist(coefs))
  }else{
    status <- 1
  }
  
  ret <- new("ModelEstimates", 
             coefficients = coefs,
             status = status,
             logL = logL
  )
  return(ret)
}

registerEstimator("spfrontier110HN",
                  new("Estimator", 
                      id = "mle", 
                      initialize = spfrontier110HN.prepare, 
                      ini.values = spfrontier110HN.ini, 
                      logL = spfrontier110HN.logL,
                      handle.estimates = spfrontier110HN.handle.estimates
                  ))