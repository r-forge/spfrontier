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
  names(ret) <- "Log-Lik SFA normal/half-normal"
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
  
  W = get("W", envir =env)
  mf <- model.frame(formula, data)
    
  y <- as.matrix(model.response(mf))
  Wy = W %*% y
  data$Wy = Wy
  
  formula = update(formula,  ~ . + Wy)
  
  sfa <- sfa.hnormal.estimator(formula, data)
  
  
  sigma <- sfa$sigma
  rm("est_success", envir=env)
  if (length(sigma)==0) {
    print("SFA is failed")
    assign("est_success", FALSE, envir=env)
    return(NULL)
  }
  
  
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

sarsfa.hnormal.estimator <- function(formula, data,W,silent=TRUE,env=new.env()){
  assign("max_ll_value", 1E10, envir=env)
  assign("W", W, envir=env)
  num = 0
  if (exists("seq_number", envir =env)){
    num = get("seq_number", envir =env)
  }
  #num = num +1 
  #assign("seq_number", num, envir=env)
  #print(paste("Run number ",num))
  estimates <- optim.estimator(formula, data, sarsfa.hnormal.lf, sarsfa.hnormal.ini, gr=sarsfa.hnormal.lf.gradient, silent=silent,env=env)
  if(!is.null(estimates)){
    est <-as.vector(estimates$estimate)
  k <- length(est)
  gamma <- as.vector(est[1:(k-3)])
  nu <- as.numeric(est[k-2])
  lambda <-  as.numeric(est[k-1])
  rho <-  as.numeric(est[k])
  sigma <- 1/nu
  beta <- gamma*sigma
  sigmaV <-sigma/sqrt(1+lambda^2)
  sigmaU <- lambda*sigmaV
  res <- list(beta=beta,sigma=sigma,lambda=lambda,rho=rho,sigmaV=sigmaV, sigmaU=sigmaU)
  return(res)
  } else{
    return(list())
  }
}