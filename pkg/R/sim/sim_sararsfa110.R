parameters <- createParDef(selection = list(n=c(100),sigmaX=10, beta0=1,beta1=-2,beta2=3, rho=0.2,rho2=0.1, sigmaV=3, sigmaU=5))

sararsfa110.dgp <- function(){
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
  u <- abs(rnorm(n, 0, sigmaU))
  y <- SpW%*%(beta0 + X %*% beta + v - u)
  dat <- data.frame(y,X)
  colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
  result <- list(formula=formula, data=dat,W=W,W2=W2)
  return(result)
}

sararsfa110.hnormal.estimator.wrapper <- function(d){
  env=new.env()
  if (exists("curenv", envir=.GlobalEnv)){
    env = get("curenv", envir=.GlobalEnv)
  }else{
    assign("curenv", env, envir=.GlobalEnv)
  }
  fake = rep(1000,7)
  
  #assign("ini_value", tv, envir=env)
  est <- sararsfa110.hnormal.estimator(d$formula,d$data,W=d$W,W2=d$W2,silent=FALSE, env)
  est_failed = exists("est_success", envir=env) && !get("est_success", envir=env)
  if (est_failed){ 
    out <- fake #Livehack
  }else{
    out <- c(est$beta,est$rho,est$rho2,est$sigmaV, est$sigmaU)
  }
  names(out) <- c("beta0_hat","beta1_hat","beta2_hat","rho","rho2", "sigmaV_hat","sigmaU_hat")
  return(out)
}

true.value <- function(){
  tv <- c(beta0, beta1, beta2,rho,rho2,sigmaV, sigmaU)
  return(tv)
}

remove("curenv",  envir=.GlobalEnv)
set.seed(0)
ezsim_sararsfa110<-ezsim(
  m             = 100,
  run           = TRUE,
  display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",rho="rho",rho2="rho2",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
  parameter_def = parameters,
  dgp           = sararsfa110.dgp,
  estimator     = sararsfa110.hnormal.estimator.wrapper,
  true_value    = true.value,
  auto_save = 100
)

runs = length(ezsim_sararsfa110$simulation_result[[1]])
for (i in 1:runs){
  if (ezsim_sararsfa110$simulation_result[[1]][[runs-i+1]][1]==1000){
    ezsim_sararsfa110$simulation_result[[1]][[runs-i+1]] = NULL
  }
}
runs = length(ezsim_sararsfa110$simulation_result[[2]])
results = data.frame()
for (i in 1:runs){
  if (ezsim_sararsfa110$simulation_result[[2]][[runs-i+1]][1]==1000){
    ezsim_sararsfa110$simulation_result[[2]][[runs-i+1]] = NULL
  }else{
    results = rbind(results,ezsim_sararsfa110$simulation_result[[2]][[runs-i+1]])
  }
}
colnames(results) <- c("beta0_hat","beta1_hat","beta2_hat","rho","rho2", "sigmaV_hat","sigmaU_hat")

ezsim_sararsfa110 = createSimulationTable(ezsim_sararsfa110)
summary(ezsim_sararsfa110)

plot(density(results$rho))
plot(density(results$rho2))
plot(density(results$mu))
# У sigmaU_hat похоже на локальный экстремум в нуле (тили не локальный, а неотличимый от отрицательной константы на таком количестве данных) - можно и их отфильтровать для презентации
plot(ezsim_sararsfa110)
plot(ezsim_sararsfa110,'density')
