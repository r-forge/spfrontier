spfrontier100HN.dgp <- function(){
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
    result <- list(formula=formula, data=dat,W=W, tv=eval(spfrontier100HN.true.value, envir=environment()))
    return(result)
}

spfrontier100HN.estimator <- function(d){
    modelEstimates <- spfrontier(d$formula,d$data,model="spfrontier100HN",W_y=d$W,logging = "info",control=list(reltol=1e-16))
    if (status(modelEstimates) > 0){ 
        fake = rep(1000,length(d$tv))
        names(fake) <- names(d$tv)
        out <- fake #Livehack
    }else{
        out <- coefficients(modelEstimates)
    }
    return(out)
}

spfrontier100HN.true.value <- function(){
    tv = c(beta0, beta1, beta2, sigmaV, sigmaU, rho)
    names(tv) <- c("Beta0","Beta1","Beta2", "SigmaV","SigmaU", "Rho")
    return(tv)
}

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
        colnames(results[[j]]) <- c("Beta0","Beta1","Beta2", "SigmaV","SigmaU", "Rho")
    }
    ezsim_ob = createSimulationTable(ezsim_ob)
    ezsim_ob$results = results
    return(ezsim_ob)
}

spfrontier100HN.test <- function(){ 
    set.seed(0)
    ezsim_spfrontier100HN<-ezsim(
        m                         = 100,
        run                     = TRUE,
        display_name    = c(beta0_hat="hat(beta)[0]",beta1_hat="hat(beta)[1]",beta2_hat="hat(beta)[2]",sigmaV_hat="hat(sigma)[v]",sigmaU_hat="hat(sigma)[u]",rho="rho"),
        parameter_def = createParDef(selection = list(n=c(100),sigmaX=10, beta0=1,beta1=-2,beta2=3,sigmaV=3, sigmaU=5, rho=0.2)),
        dgp                     = spfrontier100HN.dgp,
        estimator         = spfrontier100HN.estimator,
        true_value        = spfrontier100HN.true.value
    )
    
    ezsim_spfrontier100HN = clearFakes(ezsim_spfrontier100HN)
    
    print(summary(ezsim_spfrontier100HN))
    plot(density(ezsim_spfrontier100HN$results[[1]]$Rho))
}