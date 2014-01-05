spfrontier110HN.dgp <- function(){
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
  result <- list(formula=formula, data=dat,W=W,W2=W2, tv=eval(spfrontier110HN.true.value, envir=environment()))
  return(result)
}

spfrontier110HN.estimator <- function(d){
  modelEstimates <- spfrontier(d$formula,d$data,model="spfrontier110HN",W_y=d$W,W_v=d$W2,logging = "debug",control=list(reltol=1e-16))
  if (status(modelEstimates) > 0){ 
    fake = rep(1000,length(d$tv))
    names(fake) <- names(d$tv)
    out <- fake #Livehack
  }else{
    out <- coef(modelEstimates)
  }
  return(out)
}

spfrontier110HN.true.value <- function(){
  tv <- c(beta0, beta1, beta2,rho,rho2,sigmaV, sigmaU)
  names(tv) <- c("Beta0","Beta1","Beta2", "Rho","Rho2",  "SigmaV","SigmaU")
  return(tv)
}

ezsim_spfrontier110HN.test <- function(){ 
  set.seed(0)
  ezsim_spfrontier110HN<-ezsim(
    m             = 1,
    run           = TRUE,
    display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",rho="rho",rho2="rho2",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
    parameter_def = createParDef(selection = list(n=c(30),sigmaX=10, beta0=1,beta1=-2,beta2=3, rho=0.2,rho2=0.1, sigmaV=3, sigmaU=5)),
    dgp           = spfrontier110HN.dgp,
    estimator     = spfrontier110HN.estimator,
    true_value    = spfrontier110HN.true.value
  )
  ezsim_sararsfa110 = clearFakes(ezsim_spfrontier110HN)
  print(summary(ezsim_spfrontier110HN))
  
  plot(density(ezsim_spfrontier110HN$results[[1]]$Rho))
  plot(density(ezsim_spfrontier110HN$results[[1]]$Rho2))
}
