parameters <- createParDef(list(n=c(100),sigmaX=10, beta0=1,beta1=-2,beta2=3, sigmaV=3, sigmaU=5))

sfa.dgp <- function(){
  formula <- as.formula("y ~ X1 + X2")
  beta<-c(beta1,beta2)
  k <- length(beta)
  X <- matrix(rnorm(n*k,0,sigmaX),n, k)
  
  v <- rnorm(n,0,sigmaV)
  u <- abs(rnorm(n, 0, sigmaU))
  y <- beta0 + X %*% beta + v - u
  dat <- data.frame(y,X)
  colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
  result <- list(formula=formula, data=dat)
  return(result)
}

sfa.hnormal.estimator.wrapper <- function(d){
  est <- sfa.hnormal.estimator(d$formula,d$data,logging = "debug")
  out <- c(est$beta,est$sigmaV, est$sigmaU)
  names(out) <- c("beta0_hat","beta1_hat","beta2_hat","sigmaV_hat","sigmaU_hat")
  return(out)
}

true.value <- function(){
  tv <- c(beta0, beta1, beta2, sigmaV, sigmaU)
  return(tv)
}



set.seed(0)
ezsim_sfa<-ezsim(
  m             = 1,
  run           = TRUE,
  display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
  parameter_def = parameters,
  dgp           = sfa.dgp,
  estimator     = sfa.hnormal.estimator.wrapper,
  true_value    = true.value
)

summary(ezsim_sfa)
plot(ezsim_sfa)
plot(ezsim_sfa,'density')