sararsfa110.dgp <- function(){
  formula <- as.formula("y ~ X1 + X2")
  beta<-c(beta1,beta2)
  k <- length(beta)
  X <- matrix(rnorm(n*k,0,sigmaX),n, k)
  W <- genW(n,type="queen")
  #W = mat.or.vec(n,n)
  SpW <- solve(diag(n)-rho*W)
  print(paste("Rho2",rho2))
  W2 <- genW(n,type="rook")
  SpW2 <- solve(diag(n)-rho2*W2)
  
  mSigmaV = sigmaV*t(SpW2)%*%SpW2
  v <- rtmvnorm(1,mean = rep(0, n),sigma = mSigmaV)[1,]
  u <- abs(rnorm(n, 0, sigmaU))
  y <- SpW%*%(beta0 + X %*% beta + v - u)
  dat <- data.frame(y,X)
  colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
  result <- list(formula=formula, data=dat,W=W,W2=W2)
  return(result)
}

sararsfa110.hnormal.estimator.wrapper <- function(d){
  fake = rep(1000,7)
  est <- sararsfa110.hnormal.estimator(d$formula,d$data,W=d$W,W2=d$W2,logging="debug")
  if (!is.null(est$failed)){ 
    out <- fake #Livehack
  }else{
    out <- c(est$beta,est$rho,est$rho2,est$sigmaV, est$sigmaU)
  }
  names(out) <- c("beta0_hat","beta1_hat","beta2_hat","rho","rho2", "sigmaV_hat","sigmaU_hat")
  return(out)
}

sararsfa110.true.value <- function(){
  tv <- c(beta0, beta1, beta2,rho,rho2,sigmaV, sigmaU)
  return(tv)
}

set.seed(0)
ezsim_sararsfa110<-ezsim(
  m             = 1,
  run           = TRUE,
  display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",rho="rho",rho2="rho2",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
  parameter_def = createParDef(selection = list(n=c(30),sigmaX=10, beta0=1,beta1=-2,beta2=3, rho=0.2,rho2=0.1, sigmaV=3, sigmaU=5)),
  dgp           = sararsfa110.dgp,
  estimator     = sararsfa110.hnormal.estimator.wrapper,
  true_value    = sararsfa110.true.value
)
ezsim_sararsfa110 = clearFakes(ezsim_sararsfa110)
summary(ezsim_sararsfa110)

plot(density(ezsim_sararsfa110$results[[1]]$rho))
plot(density(ezsim_sararsfa110$results[[1]]$rho2))
plot(density(ezsim_sararsfa110$results[[1]]$mu))


# У sigmaU_hat похоже на локальный экстремум в нуле (тили не локальный, а неотличимый от отрицательной константы на таком количестве данных) - можно и их отфильтровать для презентации
plot(ezsim_sararsfa110)
plot(ezsim_sararsfa110,'density')
