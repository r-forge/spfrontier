sfa.params = function(parameters){
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

sfa.hnormal.lf <- function(parameters,env){
  y = spfrontier.env.get("y")
  X = spfrontier.env.get("X")
  p = sfa.params(parameters)
  counter = spfrontier.env.counter("sfa.hnormal.lf")
  logger.debug(paste("Evaluating 'sfa.hnormal.lf' for parameters","[run=",counter,"]:"),parameters)
  omega <- p$nu * y - X %*% p$gamma
  N <- length(y)  
  ret = Inf
  if ((p$lambda>0) && (p$nu>0)){
    ret = N * log(p$nu) - 0.5 * t(omega)%*%omega + sum(log(pnorm(-omega*p$lambda)))
  }else{
    logger.debug("'sfa.hnormal.lf' value:",ret)
  }
  if (is.nan(ret)) ret <- Inf
  logger.debug(paste("'sfa.hnormal.lf' value:",ret))
  return(-ret)
}

sfa.hnormal.lf.gradient<-function(parameters,env){
  y = spfrontier.env.get("y")
  X = spfrontier.env.get("X")
  p = sfa.params(parameters)
  logger.debug("Evaluating 'sfa.hnormal.lf.gradient' for parameters:",parameters)
  omega = p$nu * y - X %*% p$gamma
  a = -omega * p$lambda
  delta = dnorm(a)/pnorm(a)
  N = length(y)
  
  dLdGamma = t(omega)%*%X + t(delta)%*%X * p$lambda
  dLdNu = -t(omega) %*% y - t(delta)%*%y*p$lambda + N * 1/p$nu
  dLdLambda = -t(delta)%*%omega
  res = c(-dLdGamma,-dLdNu,-dLdLambda)
  names(res) = c(paste("dGamma", seq(length(dLdGamma)), sep = ""), "dNu", "dLambda")
  logger.debug("'sfa.hnormal.lf.gradient' value:",res)
  return(res)
}

sfa.hnormal.ini<-function(formula, data){
  logger.debug("Calculating initial values")
  ols <- lm(formula, data=data)
  res_ols <- resid(ols)
  sigmaU <- var(res_ols)*(1-2/pi)
  sigmaV <- var(res_ols)
  beta <- as.vector(coef(ols))
  
  par = list(beta = beta,sigmaV=sigmaV,sigmaU=sigmaU)
  p = ord.reparam(par)
  result = c(p$gamma,p$nu,p$lambda)
  logger.info("Initial values:", c(par$beta,par$sigmaV,par$sigmaU))
  return(result)
}

sfa.prepare = function(formula, data){
  mf <- model.frame(formula, data)
  y <- as.matrix(model.response(mf))
  X <- as.matrix(mf[-1])
  tm <- attr(mf, "terms")
  intercept <- attr(tm, "intercept") == 1
  if (intercept)  X <- cbind(1L,X)
  assign("X", X, envir = spfrontier.env)
  assign("y", y, envir = spfrontier.env)
}

sfa.hnormal.estimator <- function(formula, data, logging = "quiet", ini.values = NULL){
  spfrontier.env.clear()
  assign("logging.level",logging, envir = spfrontier.env)
  print(">>>>>>>>>>>>>>>>")
  logger.debug("Estimator started")
  
  sfa.prepare(formula, data)

  if (is.null(ini.values)){
    ini.values = sfa.hnormal.ini(formula, data)  
  } 
  print(ini.values)
  estimates = optim.estimator(formula, data, sfa.hnormal.lf, ini = ini.values, gr=sfa.hnormal.lf.gradient)
  res = ord.reparamBack(sfa.params(as.vector(estimates$estimate)))
  logger.info("Estimates:",c(res$beta,res$sigmaV,res$sigmaU))  
  print("")
  return(res)
}