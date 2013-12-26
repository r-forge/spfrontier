sfa.hnormal.lf <- function(parameters,env){
  y = get("y", envir =env)
  X = get("X", envir =env)
  maxValue= get("max_ll_value", envir =env)
  gamma <- parameters[1:(length(parameters)-2)]
  nu <- parameters[length(parameters)-1]
  lambda <- parameters[(length(parameters))]
  omega <- nu * y - X %*% gamma
  N <- length(y)
  if (lambda<=0) return (maxValue)
  if (nu<=0) return (maxValue)
  ret = -(N * log(nu) - 0.5 * t(omega)%*%omega + sum(log(pnorm(-omega*lambda))))
  if (is.nan(ret) || ret>maxValue) {
    ret <-maxValue
  }
  names(ret) <- "SFA normal/half-normal"
  return(ret)
}

sfa.hnormal.lf.gradient<-function(parameters,env){
  y = get("y", envir =env)
  X = get("X", envir =env)
  gamma <- parameters[1:(length(parameters)-2)]
  nu <- parameters[length(parameters)-1]
  lambda <- parameters[(length(parameters))]
  omega <- nu * y - X %*% gamma
  a <- -omega * lambda
  delta <- dnorm(a)/pnorm(a)
  N <- length(y)
  
  dLdGamma <- t(omega)%*%X + t(delta)%*%X * lambda
  dLdNu <- -t(omega) %*% y - t(delta)%*%y*lambda + N * 1/nu
  dLdLambda <- -t(delta)%*%omega
  return(c(-dLdGamma,-dLdNu,-dLdLambda))
}

sfa.hnormal.ini<-function(formula, data, env=new.env()){
  ols <- lm(formula, data=data)
  res_ols <- resid(ols)
  sigmaU <- var(res_ols)*(1-2/pi)
  sigmaV <- var(res_ols)
  
  sigma <- sqrt(sigmaV^2+sigmaU^2)
  lambda <- sigmaU/sigmaV
  nu <- 1/sigma
  beta <- as.vector(coef(ols))
  gamma <- beta/sigma
  result <-c(gamma,nu,lambda)
  names(result) <-c(paste("Gamma", seq(length(result)-2), sep = ""),"Nu","Lambda")
  return(result)
}

sfa.hnormal.estimator <- function(formula, data,silent=TRUE){
  env <- new.env()
  assign("max_ll_value", 1E10, envir=env)
  estimates <- optim.estimator(formula, data, sfa.hnormal.lf, sfa.hnormal.ini,env=env, gr=sfa.hnormal.lf.gradient, silent=silent)
  est <-as.vector(estimates$estimate)
  k <- length(est)
  gamma <- as.vector(est[1:(k-2)])
  nu <- as.numeric(est[k-1])
  lambda <-  as.numeric(est[k])
  sigma <- 1/nu
  beta <- gamma*sigma
  sigmaV <-sigma/sqrt(1+lambda^2)
  sigmaU <- lambda*sigmaV
  res <- list(beta=beta,sigma=sigma,lambda=lambda,sigmaV=sigmaV, sigmaU=sigmaU)
  return(res)
}