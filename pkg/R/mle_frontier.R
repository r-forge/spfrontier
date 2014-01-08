frontierHN.params = function(parameters){
    k = length(parameters)
    gamma = parameters[1:(k-2)]
    nu = parameters[k-1]
    lambda = parameters[k]
    names(gamma) = paste("Gamma", seq(k-2), sep = "")
    names(nu) = "Nu"
    names(lambda) = "Lambda"
    
    return(list(gamma=gamma,nu=nu,lambda=lambda))
}

olsen.reparamBack = function(ordParams){
    sigma = 1/ordParams$nu
    beta = ordParams$gamma*sigma
    sigmaV = sigma/sqrt(1+ordParams$lambda^2)
    sigmaU = ordParams$lambda*sigmaV
    
    ordParams$gamma = NULL
    ordParams$nu = NULL
    ordParams$lambda = NULL
    
    X = envirGet("X")
    names(beta) = colnames(X)
    names(sigmaV) = "SigmaV"
    names(sigmaU) = "SigmaU"
    
    return(c(list(beta=beta,sigmaV=sigmaV,sigmaU=sigmaU),ordParams))
}

olsen.reparam = function(params){
    sigma <- sqrt(params$sigmaV^2+params$sigmaU^2)
    lambda <- params$sigmaU/params$sigmaV
    nu <- 1/sigma
    gamma <- params$beta/sigma
    
    params$beta = NULL
    params$sigmaU = NULL
    params$sigmaV = NULL
    
    names(gamma) = paste("Gamma", seq(length(gamma)), sep = "")
    names(nu) = "Nu"
    names(lambda) = "Lambda"
    return(c(list(gamma=gamma,nu=nu,lambda=lambda),params))
}

frontierHN.logL <- function(parameters, quiet = F){
    if (!quiet){
        counter = envirCounter("frontierHN.logL")
        logging(paste("Evaluating 'frontierHN.logL' for parameters","[run=",counter,"]:"),parameters)
    }
    y = envirGet("y")
    X = envirGet("X")
    p = frontierHN.params(parameters)
    omega <- p$nu * y - X %*% p$gamma
    N <- length(y)    
    ret = -10e8
    if ((p$lambda>0) && (p$nu>0)){
        ret = N * log(p$nu) - 0.5 * t(omega)%*%omega + sum(log(pnorm(-omega*p$lambda)))
    }else{
        if (!quiet) logging("Parameters are out of space")
    }
    if (is.nan(ret) || (ret==-Inf)) ret = -10e8
    if (!quiet){
        logging(paste("frontierHN.logL =",-ret,"\n"))
    }
    return(-ret)
}

frontierHN.logL.gradient<-function(parameters){
    logging("Evaluating 'frontierHN.logL.gradient' for parameters:",parameters)
    y = envirGet("y")
    X = envirGet("X")
    p = frontierHN.params(parameters)
    omega = p$nu * y - X %*% p$gamma
    a = -omega * p$lambda
    delta = dnorm(a)/pnorm(a)
    N = length(y)
    
    dLdGamma = t(omega)%*%X + t(delta)%*%X * p$lambda
    dLdNu = -t(omega) %*% y - t(delta)%*%y*p$lambda + N * 1/p$nu
    dLdLambda = -t(delta)%*%omega
    res = c(-dLdGamma,-dLdNu,-dLdLambda)
    names(res) = c(paste("dGamma", seq(length(dLdGamma)), sep = ""), "dNu", "dLambda")
    return(res)
}

frontierHN.ini<-function(formula, data){
    logging("frontierHN: calculating initial values")
    ols <- lm(formula, data=data)
    res_ols <- resid(ols)
    m2 = sum(res_ols^2)/length(res_ols)
    m3 = sum(res_ols^3)/length(res_ols)
    sigmaU = (m3*sqrt(pi/2)/(1-4/pi))^(1/3)
    sigmaV = 0
    if (!is.nan(sigmaU) && (m2 - (1-2/pi)*sigmaU^2>0)){
        sigmaV = sqrt(m2 - (1-2/pi)*sigmaU^2)
    }
    if ((sigmaU<=0) || (sigmaV<=0)){
        sigmaU <- var(res_ols)*(1-2/pi)
        sigmaV <- var(res_ols)
    }
    beta <- as.vector(coef(ols))
    
    par = list(beta = beta,sigmaV=sigmaV,sigmaU=sigmaU)
    p = olsen.reparam(par)
    result = c(p$gamma,p$nu,p$lambda)
    logging("Initial values:", c(par$beta,par$sigmaV,par$sigmaU), level = "info")
    return(result)
}

prepareXY <- function(formula, data){
    mf <- model.frame(formula, data)
    y <- as.matrix(model.response(mf))
    X <- as.matrix(mf[-1])
    tm <- attr(mf, "terms")
    intercept <- attr(tm, "intercept") == 1
    if (intercept){
        X <- cbind(Intercept=1L,X)
    } 
    envirAssign("X", X)
    envirAssign("y", y)
}

frontierHN.prepare = function(formula, data,...){
    prepareXY(formula, data)
}

frontierHN.logL.direct <- function(parameters){
    k = length(parameters)
    beta = parameters[1:(k-2)]
    sigmaV = parameters[k-1]
    sigmaU = parameters[k]
    
    par <- olsen.reparam(list(beta = beta, sigmaV=sigmaV,sigmaU=sigmaU))
    pargamma <- c(par$gamma,par$nu,par$lambda)
    return(frontierHN.logL(pargamma, quiet = T))
}
frontierHN.handle.residuals <- function(estimates){
    coef <- coefficients(estimates)
    
    fittedY <- fitted(estimates)
    
    y = envirGet("y")
    resid <- y - fittedY
    rownames(resid) <- rownames(y)
    colnames(resid) <- c("Residuals")
    residuals(estimates) <- resid
    
    n <- length(resid)
    X = envirGet("X")
    
    sigma <- sqrt(coef$sigmaU^2 + coef$sigmaV^2)
    A <- resid * (coef$sigmaU / coef$sigmaV) / sigma
    u <- (dnorm(A) / (1 - pnorm(A)) - A) * coef$sigmaU * coef$sigmaV / sigma
    
    #sigmaZv <- coef$sigmaU * coef$sigmaV/sigma
    #muZv <- -resid*coef$sigmaU^2/sigma^2
    #u_h <- mtmvnorm(mean = as.vector(muZv), sigma=sigmaZv^2 * diag(n),lower=rep(0,n), doComputeVariance=FALSE)$tmean
    
    eff <- exp(-u)
    rownames(eff) <- rownames(y)
    colnames(eff) <- c("Efficiency values")
    efficiencies(estimates) <- eff
    return(estimates)
}

frontierHN.postEstimation <- function(estimates){
    logging("frontierHN: Handling estimates")
    if (status(estimates) == 0){
        coef <- olsen.reparamBack(frontierHN.params(resultParams(estimates)))
        coefficients(estimates) <- coef
        hess <- optimHess(c(coef$beta,coef$sigmaV,coef$sigmaU),frontierHN.logL.direct)
        hessian(estimates) <- hess
        
        X = envirGet("X")
        fittedY <- X %*% coef$beta
        rownames(fittedY) <- rownames(X)
        colnames(fittedY) <- c("fittedted values")
        fitted(estimates) <- fittedY
        
        estimates <- frontierHN.handle.residuals(estimates)
    }
    logging("frontierHN: Completed")
    return(estimates)
}