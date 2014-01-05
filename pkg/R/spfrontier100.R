spfrontier100HN.params = function(parameters){
  k = length(parameters)
  p = frontierHN.params(head(parameters,-1))
  p$rho = parameters[k]
  return(p)
}


spfrontier100HN.logL <- function(parameters){
  counter = envir.counter("spfrontier100HN.logL")
  logger.debug(paste("Evaluating 'spfrontier100HN.logL' for parameters","[run=",counter,"]:"),parameters)
  y = envir.get("y")
  X = envir.get("X")
  W = envir.get("W")
  p = spfrontier100HN.params(parameters)
  
  omega <- p$nu * (y-p$rho*W%*%y) - X%*%p$gamma
  N <- length(y)
  SpDet <- det(diag(N)-p$rho*W)
  ret = -10e8
  if ((p$lambda>0)&&(p$nu>0)&&(abs(p$rho)<1)&&(SpDet > 0)){
    ret = log(SpDet) + N * log(p$nu) - 0.5 * t(omega)%*%omega + sum(log(pnorm(-omega*p$lambda)))
  }else{
    logger.debug("Parameters are out of space")
  }
  if (is.nan(ret) || (ret==-Inf)) ret = -10e8
  logger.debug(paste("spfrontier100HN.logL =",-ret))
  return(-ret)
}

spfrontier100HN.logL.gradient<-function(parameters){
  logger.debug("Evaluating 'spfrontier100HN.logL.gradient' for parameters:",parameters)
  y = envir.get("y")
  X = envir.get("X")
  W = envir.get("W")
  
  p = spfrontier100HN.params(parameters)
  
  
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
  return(grad)
}

spfrontier100HN.ini <- function(formula, data){
  logger.debug("spfrontier100HN: calculating initial values")
  W = envir.get("W")
  noSpatLag = is.null(W)
  if (!noSpatLag){
    mf <- model.frame(formula, data)
    y <- as.matrix(model.response(mf))
    Wy = W %*% y
    data$Wy = Wy
    formula = update(formula,  ~ . + Wy)
  }
  sfa <- spfrontier(formula, data, model="frontierHN")
  coefs <- sfa@coefficients
  if (length(coefs$beta)==0) {
    print("SFA is failed")
    envir.assign("est.failed", TRUE)
    return(NULL)
  }
  
  beta <- coefs$beta
  if (!noSpatLag){
    coefs$rho = tail(beta, n=1)
    names(coefs$rho) = "Rho"
    coefs$beta = head(beta, -1)
  }else{
    coefs$rho = 0
  }
  p = ord.reparam(coefs)
  result = c(p$gamma,p$nu,p$lambda,p$rho)
  logger.info("Initial values:", result)
  return(result)
}

spfrontier100HN.prepare = function(formula, data,W,...){
  prepareXY(formula, data)
  envir.assign("W", W)
}

spfrontier100HN.handle.estimates <- function(estimates){
  logger.debug("spfrontier100HN: Handling estimates")
  coefs <- list()
  logL <- 0
  if (!is.null(estimates$estimate)){
    status <- 0
    coefs <- ord.reparamBack(spfrontier100HN.params(as.vector(estimates$estimate)))
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