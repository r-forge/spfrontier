parameters <- createParDef(selection = list(n=c(100),sigmaX=10, beta0=1,beta1=-2,beta2=3, rho=0.2,sigmaV=3, sigmaU=5))

sarsfa.dgp <- function(){
  formula <- as.formula("y ~ X1 + X2")
  beta<-c(beta1,beta2)
  k <- length(beta)
  X <- matrix(rnorm(n*k,0,sigmaX),n, k)
  W <- genW(n,type="queen")
  SpW <- solve(diag(n)-rho*W)
  v <- rnorm(n,0,sigmaV)
  u <- abs(rnorm(n, 0, sigmaU))
  y <- SpW%*%(beta0 + X %*% beta + v - u)
  dat <- data.frame(y,X)
  colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
  result <- list(formula=formula, data=dat,W=W)
  return(result)
}

sarsfa.hnormal.estimator.wrapper <- function(d){
  env=new.env()
  if (exists("curenv", envir=.GlobalEnv)){
    env = get("curenv", envir=.GlobalEnv)
  }else{
    assign("curenv", env, envir=.GlobalEnv)
  }
  fake = rep(1000,6)
  
  #assign("ini_value", tv, envir=env)
  est <- sarsfa.hnormal.estimator(d$formula,d$data,W=d$W,silent=TRUE, env)
  
  est_failed = exists("est_success", envir=env) && !get("est_success", envir=env)
  if (est_failed){ 
    out <- fake #Livehack
  }else{
    out <- c(est$beta,est$rho,est$sigmaV, est$sigmaU)
  }
  names(out) <- c("beta0_hat","beta1_hat","beta2_hat","rho", "sigmaV_hat","sigmaU_hat")
  return(out)
}

true.value <- function(){
  tv <- c(beta0, beta1, beta2,rho, sigmaV, sigmaU)
  return(tv)
}

remove("curenv",  envir=.GlobalEnv)
set.seed(0)
ezsim_sarsfa<-ezsim(
  m             = 100,
  run           = TRUE,
  display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",rho="rho",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
  parameter_def = parameters,
  dgp           = sarsfa.dgp,
  estimator     = sarsfa.hnormal.estimator.wrapper,
  true_value    = true.value
)

runs = length(ezsim_sarsfa$simulation_result[[1]])
for (i in 1:runs){
  if (ezsim_sarsfa$simulation_result[[1]][[runs-i+1]][1]==1000){
    ezsim_sarsfa$simulation_result[[1]][[runs-i+1]] = NULL
  }
}
runs = length(ezsim_sarsfa$simulation_result[[2]])
results = data.frame()
for (i in 1:runs){
  if (ezsim_sarsfa$simulation_result[[2]][[runs-i+1]][1]==1000){
    ezsim_sarsfa$simulation_result[[2]][[runs-i+1]] = NULL
  }else{
    results = rbind(results,ezsim_sarsfa$simulation_result[[2]][[runs-i+1]])
  }
}
colnames(results) <- c("beta0_hat","beta1_hat","beta2_hat","rho", "sigmaV_hat","sigmaU_hat")

ezsim_sarsfa = createSimulationTable(ezsim_sarsfa)
summary(ezsim_sarsfa)

plot(density(results$rho))
# У sigmaU_hat похоже на локальный экстремум в нуле (тили не локальный, а неотличимый от отрицательной константы на таком количестве данных) - можно и их отфильтровать для презентации
plot(ezsim_sarsfa)
plot(ezsim_sarsfa,'density')