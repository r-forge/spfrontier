frontierHN.params = function(parameters){
  k = length(parameters)
  gamma = parameters[1:(k-2)]
  nu = parameters[k-1]
  lambda = parameters[k]
  names(gamma) = paste("Gamma", seq(k-2), sep = "")
  names(nu) = "Nu"
  names(lambda) = "Lambda"
  
  return(list(gamma=gamma,nu=nu,lambda=lambda))
}

ord.reparamBack = function(ordParams){
  sigma = 1/ordParams$nu
  beta = ordParams$gamma*sigma
  sigmaV = sigma/sqrt(1+ordParams$lambda^2)
  sigmaU = ordParams$lambda*sigmaV
  
  ordParams$gamma = NULL
  ordParams$nu = NULL
  ordParams$lambda = NULL
  
  names(beta) = paste("Beta", seq(length(beta)), sep = "")
  names(sigmaV) = "SigmaV"
  names(sigmaU) = "SigmaU"
  
  return(c(list(beta=beta,sigmaV=sigmaV,sigmaU=sigmaU),ordParams))
}

ord.reparam = function(params){
  sigma <- sqrt(params$sigmaV^2+params$sigmaU^2)
  lambda <- params$sigmaU/params$sigmaV
  nu <- 1/sigma
  gamma <- params$beta/sigma
  
  params$beta = NULL
  params$sigmaU = NULL
  params$sigmaV = NULL
  
  names(gamma) = paste("Gamma", seq(length(gamma)), sep = "")
  names(nu) = "Nu"
  names(lambda) = "Lambda"
  return(c(list(gamma=gamma,nu=nu,lambda=lambda),params))
}

frontierHN.logL <- function(parameters,env){
  counter = envir.counter("frontierHN.logL")
  logger.debug(paste("Evaluating 'frontierHN.logL' for parameters","[run=",counter,"]:"),parameters)
  y = envir.get("y")
  X = envir.get("X")
  p = frontierHN.params(parameters)
  omega <- p$nu * y - X %*% p$gamma
  N <- length(y)  
  ret = -10e8
  if ((p$lambda>0) && (p$nu>0)){
    ret = N * log(p$nu) - 0.5 * t(omega)%*%omega + sum(log(pnorm(-omega*p$lambda)))
  }else{
    logger.debug("Parameters are out of space")
  }
  if (is.nan(ret) || (ret==-Inf)) ret = -10e8
  logger.debug(paste("frontierHN.logL = ",-ret,sep=""))
  return(-ret)
}

frontierHN.logL.gradient<-function(parameters,env){
  logger.debug("Evaluating 'frontierHN.logL.gradient' for parameters:",parameters)
  y = envir.get("y")
  X = envir.get("X")
  p = frontierHN.params(parameters)
  omega = p$nu * y - X %*% p$gamma
  a = -omega * p$lambda
  delta = dnorm(a)/pnorm(a)
  N = length(y)
  
  dLdGamma = t(omega)%*%X + t(delta)%*%X * p$lambda
  dLdNu = -t(omega) %*% y - t(delta)%*%y*p$lambda + N * 1/p$nu
  dLdLambda = -t(delta)%*%omega
  res = c(-dLdGamma,-dLdNu,-dLdLambda)
  names(res) = c(paste("dGamma", seq(length(dLdGamma)), sep = ""), "dNu", "dLambda")
  return(res)
}

frontierHN.ini<-function(formula, data){
  logger.debug("frontierHN: calculating initial values")
  ols <- lm(formula, data=data)
  res_ols <- resid(ols)
  m2 = sum(res_ols^2)/length(res_ols)
  m3 = sum(res_ols^3)/length(res_ols)
  sigmaU = (m3*sqrt(pi/2)/(1-4/pi))^(1/3)
  sigmaV = 0
  if (!is.nan(sigmaU) && (m2 - (1-2/pi)*sigmaU^2>0)){
    sigmaV = sqrt(m2 - (1-2/pi)*sigmaU^2)
  }
  if ((sigmaU<=0) || (sigmaV<=0)){
    sigmaU <- var(res_ols)*(1-2/pi)
    sigmaV <- var(res_ols)
  }
  beta <- as.vector(coef(ols))
  
  par = list(beta = beta,sigmaV=sigmaV,sigmaU=sigmaU)
  p = ord.reparam(par)
  result = c(p$gamma,p$nu,p$lambda)
  logger.info("Initial values:", c(par$beta,par$sigmaV,par$sigmaU))
  return(result)
}

prepareXY <- function(formula, data){
  mf <- model.frame(formula, data)
  y <- as.matrix(model.response(mf))
  X <- as.matrix(mf[-1])
  tm <- attr(mf, "terms")
  intercept <- attr(tm, "intercept") == 1
  if (intercept)  X <- cbind(1L,X)
  envir.assign("X", X)
  envir.assign("y", y)
}

frontierHN.prepare = function(formula, data,...){
  prepareXY(formula, data)
}

frontierHN.handle.estimates <- function(estimates){
  logger.debug("frontierHN: Handling estimates")
  coefs <- list()
  logL <- 0
  if (!is.null(estimates$estimate)){
    status <- 0
    coefs <- ord.reparamBack(frontierHN.params(as.vector(estimates$estimate)))
    logL <- estimates$value
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

registerEstimator("frontierHN",
                  new("Estimator", 
                      id = "mle", 
                      initialize = frontierHN.prepare, 
                      ini.values = frontierHN.ini, 
                      logL = frontierHN.logL, 
                      gradient = frontierHN.logL.gradient, 
                      handle.estimates = frontierHN.handle.estimates
                      ))