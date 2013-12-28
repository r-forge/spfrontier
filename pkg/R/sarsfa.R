sarsfa.params = function(parameters){
  k = length(parameters)
  p = sfa.params(head(parameters,-1))
  p$Rho = parameters[k]
  p$Rho = rho
  return(p)
}


sarsfa.hnormal.lf <- function(parameters,env){
  y = get("y", envir =env)
  X = get("X", envir =env)
  W = get("W", envir =env)
  maxValue= get("max_ll_value", envir =env)
  k = length(parameters)
  
  gamma <- parameters[1:(k-3)]
  nu <- parameters[k-2]
  lambda <- parameters[k-1]
  
  
  rho <- parameters[k]
  
  omega <- nu * (y-rho*W%*%y) - X%*%gamma
  N <- length(y)
  SpDet <- log(det(diag(N)-rho*W))
  if (lambda<=0) return (maxValue)
  if (nu<=0) return (maxValue)
  if (abs(rho)>1) return (maxValue)
  ret = -(SpDet + N * log(nu) - 0.5 * t(omega)%*%omega + sum(log(pnorm(-omega*lambda))))
  if (is.nan(ret) || ret>maxValue) {
    ret <-maxValue
  }
  return(ret)
}

sarsfa.hnormal.lf.gradient<-function(parameters,env){
  
  y = get("y", envir =env)
  X = get("X", envir =env)
  W = get("W", envir =env)
  k = length(parameters)
  gamma <- parameters[1:(k-3)]
  nu <- parameters[k-2]
  lambda <- parameters[k-1]
  rho <- parameters[k]
  
  omega <- nu * (y-rho*W%*%y) - X%*%gamma
  a <- -omega * lambda
  delta <- dnorm(a)/pnorm(a)
  N <- length(y)
  
  dLdGamma <- t(omega)%*%X + t(delta)%*%X * lambda
  dLdNu <- -t(omega)%*%y - t(delta)%*%y*lambda + N * 1/nu
  dLdLambda <- -t(delta)%*%omega
  mat = solve(diag(N)-rho*W) %*% (-W)
  dLdRho = nu*(t(omega)+lambda*t(delta))%*%W%*%y+sum(diag(mat))
  #dLdRho = abs(rho-0.2)
  grad = c(-dLdGamma,-dLdNu,-dLdLambda,-dLdRho)
  return(grad)
}

sarsfa.hnormal.ini<-function(formula, data, env=new.env()){
  logger.debug("Calculating initial values")
  
  W = spfrontier.env.get("W")
  mf <- model.frame(formula, data)
  y <- as.matrix(model.response(mf))
  Wy = W %*% y
  data$Wy = Wy
  formula = update(formula,  ~ . + Wy)
  
  sfa <- sfa.hnormal.estimator(formula, data)
  
  
  sigmaV <- sfa$sigmaV
  if (length(sigmaV)==0) {
    print("SFA is failed")
    assign("est_success", FALSE, envir=env)
    return(NULL)
  }
  
  
  beta <- sfa$beta
  rho = tail(beta, n=1)
  sfa$beta = head(beta, -1)
  
  p = ord.reparam(sfa)
  result = c(p$gamma,p$nu,p$lambda,p$rho)
  logger.info("Initial values:", result)
  return(result)
  
  lambda <- sfa$lambda
  nu <- 1/sigma
  
  beta <- sfa$beta
  rho = tail(beta, n=1)
  gamma <- beta/sigma
  gamma = head(gamma, -1)
  result <-c(gamma,nu,lambda,rho)
  names(result) <-c(paste("Gamma", seq(length(gamma)), sep = ""),"Nu","Lambda", "Rho")
  
  return(result)
}

sarsfa.hnormal.estimator <- function(formula, data,logging = "quiet", ini.values = NULL){
  assign("logging.level",logging, envir = spfrontier.env)
  print(">>>>>>>>>>>>>>>>")
  logger.debug("Estimator started")
  sfa.prepare(formula, data)
  if (is.null(ini.values)){
    ini.values = sarsfa.hnormal.ini(formula, data)  
  } 
  
  assign("W", W, envir=spfrontier.env)
  
  estimates = optim.estimator(formula, data, sarsfa.hnormal.lf, ini.values, gr=sarsfa.hnormal.lf.gradient)
  res = ord.reparamBack(sarsfa.params(as.vector(estimates$estimate)))
  logger.info("Estimates:",c(res$beta,res$sigmaV,res$sigmaU,res$rho))  
  print("")
  return(res)
}