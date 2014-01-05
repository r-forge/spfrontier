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
  fake = rep(1000,6)
  est <- sarsfa.hnormal.estimator(d$formula,d$data,W=d$W,logging="info")
  if (!is.null(est$failed)){ 
    out <- fake #Livehack
  }else{
    out <- c(est$beta,est$rho,est$sigmaV, est$sigmaU)
  }
  names(out) <- c("beta0_hat","beta1_hat","beta2_hat","rho", "sigmaV_hat","sigmaU_hat")
  return(out)
}

sarsfa.true.value <- function(){
  tv <- c(beta0, beta1, beta2,rho, sigmaV, sigmaU)
  return(tv)
}

set.seed(0)
ezsim_sarsfa<-ezsim(
  m             = 100,
  run           = TRUE,
  display_name  = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",rho="rho",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]"),
  parameter_def = createParDef(selection = list(n=c(100),sigmaX=10, beta0=1,beta1=-2,beta2=3, rho=0.2,sigmaV=3, sigmaU=5)),
  dgp           = sarsfa.dgp,
  estimator     = sarsfa.hnormal.estimator.wrapper,
  true_value    = sarsfa.true.value
)

clearFakes = function(ezsim_ob){
  parSets = length(ezsim_ob$simulation_result)
  results = list()
  for (j in 1:parSets){
    results[[j]] = data.frame()
    runs = length(ezsim_ob$simulation_result[[j]])
    for (i in 1:runs){
      if (ezsim_ob$simulation_result[[j]][[runs-i+1]][1]==1000){
        ezsim_ob$simulation_result[[j]][[runs-i+1]] = NULL
      }else{
        results[[j]] = rbind(results[[j]],ezsim_ob$simulation_result[[j]][[runs-i+1]])
      }
    }
    colnames(results[[j]]) <- c("beta0_hat","beta1_hat","beta2_hat","rho", "sigmaV_hat","sigmaU_hat")
  }
  ezsim_ob = createSimulationTable(ezsim_ob)
  ezsim_ob$results = results
  return(ezsim_ob)
}

ezsim_sarsfa = clearFakes(ezsim_sarsfa)
summary(ezsim_sarsfa)

plot(density(ezsim_sarsfa$results[[1]]$rho))
# У sigmaU_hat похоже на локальный экстремум в нуле (тили не локальный, а неотличимый от отрицательной константы на таком количестве данных) - можно и их отфильтровать для презентации
plot(ezsim_sarsfa)
plot(ezsim_sarsfa,'density')