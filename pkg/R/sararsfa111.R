sararsfa111.params = function(parameters){
  k = length(parameters)
  p = sararsfa110.params(head(parameters,-2))
  p$rho3 = parameters[k-1]
  p$mu = parameters[k]
  return(p)
}

sararsfa111.tnormal.lf <- function(parameters){
  y = envir.get("y")
  X = envir.get("X")
  W = envir.get("W")
  W2 = envir.get("W2")
  W3 = envir.get("W3")
  p = sararsfa111.params(parameters)
  counter = envir.counter("sararsfa111.hnormal.lf")
  logger.debug(paste("Evaluating 'sararsfa111.tnormal.lf' for parameters","[run=",counter,"]:"),parameters) 
  ret = -10e8
  
  
  if ((p$sigmaV>0) && (p$sigmaU>0) && (abs(p$rho)<1) && (abs(p$rho2)<1) && (abs(p$rho3)<1)){
      n = length(y)
      e <- y  - p$rho * W %*% y -  X %*% p$beta
      I <- diag(n)
      
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
        ret<- log(det(I-p$rho*W))-log(ptmvnorm(lowerx=rep(-Inf, n),upperx = rep(0, n),mean=as.vector(t(vMu)), sigma=mOmega))+log(dmvnorm(x=as.vector(e+vMu),mean=rep(0, n), sigma=mC))+log(ptmvnorm(lowerx=rep(-Inf, n),upperx = rep(0, n),mean=as.vector(t(-mA%*%e-mD%*%vMu)), sigma=mB))
      }, error = function(e){
        logger.warn(e$message)
      })
  }else{
    logger.debug("Parameters are out of space")
  }
  if (is.nan(ret) || (ret==-Inf)) ret = -10e8
  logger.debug(paste("'sararsfa111.hnormal.lf' value:",-ret))
  return(-ret)
}

sararsfa111.hnormal.lf <- function(parameters){
  return(sararsfa111.tnormal.lf(c(parameters,0)))
}
sararsfa111.hnormal.ini<-function(formula, data){
  ini = sararsfa111.tnormal.ini(formula, data)
  return(head(ini, -1))
}

sararsfa111.tnormal.ini<-function(formula, data){
  logger.debug("Calculating initial values")
  W = envir.get("W")
  sarsfa <- sarsfa.hnormal.estimator(formula, data, W)
  if (length(sarsfa$beta)==0) {
    print("sarsfa is failed")
    envir.assign("est.failed", TRUE)
    return(NULL)
  }
  
  y = envir.get("y")
  X = envir.get("X")
  W2 = envir.get("W2")
  e = y - sarsfa$rho * y - X %*% sarsfa$beta
  W2e = W2 %*% e
  model = lm(e ~ W2e-1, data = data.frame(e, W2e))

  
  rho2 = coef(model)
  rho3 = 0
  ve = resid(model)
  mu = - mean(ve)
  result <-c(sarsfa$beta,sarsfa$sigmaV,sarsfa$sigmaU,sarsfa$rho, rho2, rho3, mu)
  names(result) <-c(paste("Beta", seq(from=0, to=length(sarsfa$beta)-1),sep = ""),"SigmaV","SigmaU", "Rho", "Rho2", "Rho3", "Mu")
  logger.info("Initial values:", result)
  return(result)
}

sararsfa111.tnormal.estimator <- function(formula, data,W,W2,W3,logging = "quiet", ini.values = NULL,control=NULL){
  envir.init()
  envir.assign("logging.level",logging)
  if (!is.null(control)) envir.assign("control",control)
  logger.start()
  logger.debug("Estimator started")
  
  sfa.prepare(formula, data)
  envir.assign("W", W)
  envir.assign("W2", W2)
  envir.assign("W3", W3)
  
  if (is.null(ini.values)){
    ini.values = sararsfa111.hnormal.ini(formula, data)  
  } 
  
  estimates <- optim.estimator(formula, data, sararsfa111.hnormal.lf, ini.values, gr=NULL)
  est_failed = !is.null(envir.get("est.failed"))
  res = rep(1000,length(ini.values))
  if(est_failed){
    res$failed = T
    logger.info("Estimator failed")  
  }else{
    res = sararsfa111.params(estimates$estimate)
    logger.info("Estimates:",c(res$beta,res$sigmaV,res$sigmaU,res$rho,res$rho2,res$rho3,res$mu))  
  }
  print("")
  envir.finalise()
  return(res)
}