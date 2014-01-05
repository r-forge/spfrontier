frontierHN.dgp <- function(){
  formula <- as.formula("y ~ X1 + X2")
  beta<-c(beta1,beta2)
  
  k <- length(beta)
  X <- matrix(rnorm(n*k,0,sigmaX),n, k)
  
  v <- rnorm(n,0,sigmaV)
  u <- abs(rnorm(n, 0, sigmaU))
  y <- beta0 + X %*% beta + v - u
  dat <- data.frame(y,X)
  colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
  result <- list(formula=formula, data=dat, tv=eval(frontierHN.true.value, envir=environment()))
  return(result)
}

frontierHN.estimator <- function(d){
  modelEstimates <- spfrontier(d$formula,d$data,model="frontierHN", logging = "debug")
  out <- coef(modelEstimates)
  return(out)
}

frontierHN.true.value <- function(){
  tv = c(beta0, beta1, beta2, sigmaV, sigmaU)
  names(tv) <- c("Beta0","Beta1","Beta2", "SigmaV","SigmaU")
  return(tv)
}

frontierHN.test <- function(){
  set.seed(0)
  ezsim_frontierHN <- ezsim(
    m             = 10,
    run           = TRUE,
    display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
    parameter_def = createParDef(list(n=c(100),sigmaX=10, beta0=1,beta1=-2,beta2=3, sigmaV=3, sigmaU=5)),
    dgp           = frontierHN.dgp,
    estimator     = frontierHN.estimator,
    true_value    = frontierHN.true.value
  )
  
  print(summary(ezsim_frontierHN))
  #plot(ezsim_frontierHN)
  #plot(ezsim_frontierHN,'density')
}