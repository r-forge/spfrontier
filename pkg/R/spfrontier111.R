spfrontier111TN.params = function(parameters){
  k = length(parameters)
  p = spfrontier110HN.params(head(parameters,-2))
  p$rho3 = parameters[k-1]
  p$mu = parameters[k]
  return(p)
}

spfrontier111TN.logL <- function(parameters){
  counter = envir.counter("spfrontier111TN.logL")
  logger.debug(paste("Evaluating 'spfrontier111TN.logL' for parameters","[run=",counter,"]:"),parameters) 
  y = envir.get("y")
  X = envir.get("X")
  W = envir.get("W")
  W2 = envir.get("W2")
  W3 = envir.get("W3")
  p = spfrontier111TN.params(parameters)
  ret = -10e8
  n = length(y)
  I <- diag(n)
  SpDet <- det(I-p$rho*W)
  
  if ((p$sigmaV>0) && (p$sigmaU>0) && (abs(p$rho)<1) && (abs(p$rho2)<1) && (abs(p$rho3)<1) && (SpDet>0)){
      
      e <- y  - p$rho * W %*% y -  X %*% p$beta
      
      Sp2 <- solve(I-p$rho2*W2)
      Sp3 <- solve(I-p$rho3*W3)
      
      mSigma = p$sigmaV^2*t(Sp2)%*%Sp2
      mOmega = p$sigmaU^2*t(Sp3)%*%Sp3
      mC = mSigma + mOmega
      mB = mOmega%*%solve(mC)%*%mSigma 
      mA = mB %*% solve(mSigma)
      mD = -mB %*% solve(mOmega)
      vMu = rep(p$mu, n)
      tryCatch({
        ret<- log(SpDet)-log(ptmvnorm(lowerx=rep(-Inf, n),upperx = rep(0, n),mean=as.vector(t(vMu)), sigma=mOmega))+log(dmvnorm(x=as.vector(e+vMu),mean=rep(0, n), sigma=mC))+log(ptmvnorm(lowerx=rep(-Inf, n),upperx = rep(0, n),mean=as.vector(t(-mA%*%e-mD%*%vMu)), sigma=mB))
      }, error = function(e){
        logger.warn(e$message)
      })
  }else{
    logger.debug("Parameters are out of space")
  }
  if (is.nan(ret) || (ret==-Inf)) ret = -10e8
  logger.debug(paste("spfrontier111TN.logL =",-ret))
  return(-ret)
}

#sararsfa111.hnormal.lf <- function(parameters){
#  return(sararsfa111.tnormal.lf(c(parameters,0)))
#}
#sararsfa111.hnormal.ini<-function(formula, data){
#  ini = sararsfa111.tnormal.ini(formula, data)
#  return(head(ini, -1))
#}

spfrontier111TN.ini<-function(formula, data){
  logger.debug("spfrontier111TN: calculating initial values")
  W = envir.get("W")
  iniModel <- spfrontier(formula,data,model="spfrontier100HN",W_y=W)
  coefs <- iniModel@coefficients
  if (status(iniModel) > 0){
    iniModel <- spfrontier(formula,data,model="frontierHN")
    coefs <- iniModel@coefficients
    coefs$rho <- 0
  }
  
  y = envir.get("y")
  X = envir.get("X")
  W2 = envir.get("W2")
  e = y - coefs$rho * y - X %*% coefs$beta
  W2e = W2 %*% e
  model = lm(e ~ W2e-1, data = data.frame(e, W2e))

  
  rho2 = coef(model)
  rho3 = 0
  ve = resid(model)
  mu = - mean(ve)
  result <-c(coefs$beta,coefs$sigmaV,coefs$sigmaU,coefs$rho, rho2, rho3, mu)
  names(result) <-c(paste("Beta", seq(from=0, to=length(coefs$beta)-1),sep = ""),"SigmaV","SigmaU", "Rho", "Rho2", "Rho3", "Mu")
  logger.info("Initial values:", result)
  return(result)
}

spfrontier111TN.prepare = function(formula, data,W,W2,W3,...){
  prepareXY(formula, data)
  envir.assign("W", W)
  envir.assign("W2", W2)
  envir.assign("W3", W3)
}

spfrontier111TN.handle.estimates <- function(estimates){
  logger.debug("spfrontier111TN: Handling estimates")
  coefs <- list()
  logL <- 0
  if (!is.null(estimates$estimate)){
    status <- 0
    coefs <- spfrontier111TN.params(as.vector(estimates$estimate))
    
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

registerEstimator("spfrontier111TN",
                  new("Estimator", 
                      id = "mle", 
                      initialize = spfrontier111TN.prepare, 
                      ini.values = spfrontier111TN.ini, 
                      logL = spfrontier111TN.logL,
                      handle.estimates = spfrontier111TN.handle.estimates
                  ))