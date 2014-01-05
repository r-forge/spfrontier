sararsfa111.dgp <- function(){
  print("sararsfa111.dgp")
  formula <- as.formula("y ~ X1 + X2")
  beta<-c(beta1,beta2)
  k <- length(beta)
  X <- matrix(rnorm(n*k,0,sigmaX),n, k)
  W <- genW(n,type="queen")
  SpW <- solve(diag(n)-rho*W)
  
  W2 <- genW(n,type="rook")
  SpW2 <- solve(diag(n)-rho2*W2)
  mSigmaV = sigmaV*t(SpW2)%*%SpW2
  v <- rtmvnorm(1,mean = rep(0, n),sigma = mSigmaV)[1,]
  
  
  W3 <- genW(n,type="queen")
  SpW3 <- solve(diag(n)-rho3*W3)
  mSigmaU = sigmaU*t(SpW3)%*%SpW3
  
  #Вот так долго - но, возможно, правильно
  #u <- rtmvnorm2(1,mean = rep(mu, n),sigma = mSigmaU, lower=rep(0, n))[1,]
  
  print(mu)
  u <- abs(rtmvnorm(1,mean = rep(mu, n),sigma = mSigmaU)[1,])
  
  y <- SpW%*%(beta0 + X %*% beta + v - u)
  dat <- data.frame(y,X)
  colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
  
  tv = c(beta0, beta1, beta2, sigmaV, sigmaU, rho, rho2, rho3, mu)
  names(tv) <- c("Beta0","Beta1","Beta2", "SigmaV","SigmaU","rho","rho2","rho3","mu")
  
  result <- list(formula=formula, data=dat,W=W,W2=W2,W3=W3,tv=tv)
  
  print("DGP generated")
  return(result)
}

sararsfa111.tnormal.estimator.wrapper <- function(d){
  print(paste("logLik at true value",sararsfa111.tnormal.lf(d$tv)))
  fake = rep(1000,9)
  est <- sararsfa111.tnormal.estimator(d$formula,d$data,W=d$W,W2=d$W2,W3=d$W3,logging="debug")
  if (!is.null(est$failed)){  
    out <- fake #Livehack
  }else{
    out <- c(est$beta,est$rho,est$rho2,est$rho3,est$mu,est$sigmaV, est$sigmaU)
  }
  names(out) <- c("beta0_hat","beta1_hat","beta2_hat","rho","rho2","rho3","mu", "sigmaV_hat","sigmaU_hat")
  return(out)
}

sararsfa111.true.value <- function(){
  tv <- c(beta0, beta1, beta2,rho,rho2,rho3,mu,sigmaV, sigmaU)
  return(tv)
}

set.seed(0)
ezsim_sararsfa111<-ezsim(
  m             = 1,
  run           = TRUE,
  display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",rho="rho",rho2="rho2",rho3="rho3",mu="mu",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
  parameter_def = createParDef(selection = list(n=c(300),sigmaX=10, beta0=1,beta1=-2,beta2=3, rho=0.3,rho2=0.2,rho3=0.1, mu=1, sigmaV=3, sigmaU=5)),
  dgp           = sararsfa111.dgp,
  estimator     = sararsfa111.tnormal.estimator.wrapper,
  true_value    = sararsfa111.true.value
)



ezsim_sararsfa111 <- clearFakes(ezsim_sararsfa111)
summary(ezsim_sararsfa111)

plot(density(ezsim_sararsfa111$results[[1]]$rho))
plot(density(ezsim_sararsfa111$results[[1]]$rho2))
plot(density(ezsim_sararsfa111$results[[1]]$mu))

plot(ezsim_sararsfa111)
plot(ezsim_sararsfa111,'density')