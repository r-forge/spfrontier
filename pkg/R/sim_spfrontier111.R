spfrontier111TN.dgp <- function(){
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
  
  u <- abs(rtmvnorm(1,mean = rep(mu, n),sigma = mSigmaU)[1,])
  
  y <- SpW%*%(beta0 + X %*% beta + v - u)
  dat <- data.frame(y,X)
  colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
  
  
  result <- list(formula=formula, data=dat,W=W,W2=W2,W3=W3, tv=eval(spfrontier110HN.true.value, envir=environment()))
  
  print("DGP generated")
  return(result)
}

spfrontier111TN.estimator <- function(d){
  #print(paste("logLik at true value",sararsfa111.tnormal.lf(d$tv)))
  modelEstimates <- spfrontier(d$formula,d$data,model="spfrontier111TN",W_y=d$W,W_v=d$W2,W_u=d$W3,logging = "debug",control=list(reltol=1e-16))
  if (status(modelEstimates) > 0){ 
    fake = rep(1000,length(d$tv))
    names(fake) <- names(d$tv)
    out <- fake #Livehack
  }else{
    out <- coef(modelEstimates)
  }
  return(out)
}

spfrontier111TN.true.value <- function(){
  tv = c(beta0, beta1, beta2, rho, rho2, rho3, mu, sigmaV, sigmaU)
  names(tv) <- c("Beta0","Beta1","Beta2","rho","rho2","rho3","mu", "SigmaV","SigmaU")
  return(tv)
}

spfrontier111TN.test <- function(){
  set.seed(0)
  ezsim_spfrontier111TN<-ezsim(
    m             = 1,
    run           = TRUE,
    display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",rho="rho",rho2="rho2",rho3="rho3",mu="mu",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
    parameter_def = createParDef(selection = list(n=c(300),sigmaX=10, beta0=1,beta1=-2,beta2=3, rho=0.3,rho2=0.2,rho3=0.1, mu=1, sigmaV=3, sigmaU=5)),
    dgp           = spfrontier111TN.dgp,
    estimator     = spfrontier111TN.estimator,
    true_value    = spfrontier111TN.true.value
  )
  
  
  
  ezsim_spfrontier111TN <- clearFakes(ezsim_spfrontier111TN)
  summary(ezsim_sararsfa111)
  
  plot(density(ezsim_spfrontier111TN$results[[1]]$rho))
  plot(density(ezsim_spfrontier111TN$results[[1]]$rho2))
  plot(density(ezsim_spfrontier111TN$results[[1]]$ru))
}