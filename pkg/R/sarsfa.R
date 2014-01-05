sarsfa.params = function(parameters){
  k = length(parameters)
  p = sfa.params(head(parameters,-1))
  p$rho = parameters[k]
  return(p)
}


sarsfa.hnormal.lf <- function(parameters){
  y = envir.get("y")
  X = envir.get("X")
  W = envir.get("W")
  p = sarsfa.params(parameters)
  counter = envir.counter("sarsfa.hnormal.lf")
  logger.debug(paste("Evaluating 'sarsfa.hnormal.lf' for parameters","[run=",counter,"]:"),parameters)
  
  omega <- p$nu * (y-p$rho*W%*%y) - X%*%p$gamma
  N <- length(y)
  SpDet <- log(det(diag(N)-p$rho*W))
  ret = -10e8
  if ((p$lambda>0)&&(p$nu>0)&&(abs(p$rho)<1)){
    ret = SpDet + N * log(p$nu) - 0.5 * t(omega)%*%omega + sum(log(pnorm(-omega*p$lambda)))
  }else{
    logger.debug("Parameters are out of space")
  }
  if (is.nan(ret) || (ret==-Inf)) ret = -10e8
  logger.debug(paste("'sarsfa.hnormal.lf' value:",-ret))
  return(-ret)
}

sarsfa.hnormal.lf.gradient<-function(parameters){
  y = envir.get("y")
  X = envir.get("X")
  W = envir.get("W")
  
  p = sarsfa.params(parameters)
  logger.debug("Evaluating 'sarsfa.hnormal.lf.gradient' for parameters:",parameters)
  
  
  omega <- p$nu * (y-p$rho*W%*%y) - X%*%p$gamma
  a <- -omega * p$lambda
  delta <- dnorm(a)/pnorm(a)
  N <- length(y)
  
  dLdGamma <- t(omega)%*%X + t(delta)%*%X * p$lambda
  dLdNu <- -t(omega)%*%y - t(delta)%*%y*p$lambda + N * 1/p$nu
  dLdLambda <- -t(delta)%*%omega
  mat = solve(diag(N)-p$rho*W) %*% (-W)
  dLdRho = p$nu*(t(omega)+p$lambda*t(delta))%*%W%*%y+sum(diag(mat))
  #dLdRho = abs(rho-0.2)
  grad = c(-dLdGamma,-dLdNu,-dLdLambda,-dLdRho)
  names(grad) = c(paste("dGamma", seq(length(dLdGamma)), sep = ""), "dNu", "dLambda", "dRho")
  logger.debug("'sarsfa.hnormal.lf.gradient' value:",grad)
  return(grad)
}

sarsfa.hnormal.ini<-function(formula, data){
  logger.debug("Calculating initial values")
  W = envir.get("W")
  noSpatLag = all(mat.or.vec(10,10) == 0)
  if (!noSpatLag){
    mf <- model.frame(formula, data)
    y <- as.matrix(model.response(mf))
    Wy = W %*% y
    data$Wy = Wy
    formula = update(formula,  ~ . + Wy)
  }
  sfa <- sfa.hnormal.estimator(formula, data)
  
  if (length(sfa$beta)==0) {
    print("SFA is failed")
    envir.assign("est.failed", TRUE)
    return(NULL)
  }
  
  beta <- sfa$beta
  if (!noSpatLag){
    sfa$rho = tail(beta, n=1)
    names(sfa$rho) = "Rho"
    sfa$beta = head(beta, -1)
  }else{
    sfa$rho = 0
  }
  p = ord.reparam(sfa)
  result = c(p$gamma,p$nu,p$lambda,p$rho)
  logger.info("Initial values:", result)
  return(result)
}

sarsfa.hnormal.estimator <- function(formula, data,W,logging = "quiet", ini.values = NULL,control=NULL){
  envir.init()
  envir.assign("logging.level",logging)
  if (!is.null(control)) envir.assign("control",control)
  logger.start()
  logger.debug("Estimator started")
  
  sfa.prepare(formula, data)
  envir.assign("W", W)
  
  if (is.null(ini.values)){
    ini.values = sarsfa.hnormal.ini(formula, data)  
  } 
  
  estimates = optim.estimator(formula, data, sarsfa.hnormal.lf, ini.values, gr=sarsfa.hnormal.lf.gradient)
  est_failed = !is.null(envir.get("est.failed"))
  res = ord.reparamBack(sarsfa.params(estimates$estimate))
  if(est_failed){
    res$failed = T
  }
  logger.info("Estimates:",c(res$beta,res$sigmaV,res$sigmaU,res$rho))  
  print("")
  envir.finalise()
  return(res)
}