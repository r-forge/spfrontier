spfrontier100HN.params = function(parameters){
    k = length(parameters)
    p = frontierHN.params(head(parameters,-1))
    p$rho = parameters[k]
    return(p)
}


spfrontier100HN.logL <- function(parameters, quiet = F){
    if (!quiet){
        counter = envirCounter("spfrontier100HN.logL")
        logging(paste("Evaluating 'spfrontier100HN.logL' for parameters","[run=",counter,"]:"),parameters)
    }
    y = envirGet("y")
    X = envirGet("X")
    W = envirGet("W")
    p = spfrontier100HN.params(parameters)
    
    omega <- p$nu * (y-p$rho*W%*%y) - X%*%p$gamma
    N <- length(y)
    SpDet <- det(diag(N)-p$rho*W)
    ret = -10e8
    if ((p$lambda>0)&&(p$nu>0)&&(abs(p$rho)<1)&&(SpDet > 0)){
        ret = log(SpDet) + N * log(p$nu) - 0.5 * t(omega)%*%omega + sum(log(pnorm(-omega*p$lambda)))
    }else{
        if (!quiet) logging("Parameters are out of space")
    }
    if (is.nan(ret) || (ret==-Inf)) ret = -10e8
    if (!quiet){
        logging(paste("spfrontier100HN.logL =",-ret,"\n"))
    }
    return(-ret)
}

spfrontier100HN.logL.gradient<-function(parameters){
    logging("Evaluating 'spfrontier100HN.logL.gradient' for parameters:",parameters)
    y = envirGet("y")
    X = envirGet("X")
    W = envirGet("W")
    
    p = spfrontier100HN.params(parameters)
    
    
    omega <- p$nu * (y-p$rho*W%*%y) - X%*%p$gamma
    a <- -omega * p$lambda
    delta <- dnorm(a)/pnorm(a)
    N <- length(y)
    
    dLdGamma <- t(omega)%*%X + t(delta)%*%X * p$lambda
    dLdNu <- -t(omega)%*%y - t(delta)%*%y*p$lambda + N * 1/p$nu
    dLdLambda <- -t(delta)%*%omega
    mat = solve(diag(N)-p$rho*W) %*% (-W)
    dLdRho = p$nu*(t(omega)+p$lambda*t(delta))%*%W%*%y+sum(diag(mat))
    #dLdRho = abs(rho-0.2)
    grad = c(-dLdGamma,-dLdNu,-dLdLambda,-dLdRho)
    names(grad) = c(paste("dGamma", seq(length(dLdGamma)), sep = ""), "dNu", "dLambda", "dRho")
    return(grad)
}

spfrontier100HN.ini <- function(formula, data){
    logging("spfrontier100HN: calculating initial values")
    W = envirGet("W")
    noSpatLag = is.null(W)
    if (!noSpatLag){
        mf <- model.frame(formula, data)
        y <- as.matrix(model.response(mf))
        Wy = W %*% y
        data$Wy = Wy
        formula = update(formula,    ~ . + Wy)
    }
    sfa <- spfrontier(formula, data, model="frontierHN")
    coefs <- sfa@coefficients
    if (length(coefs$beta)==0) {
        print("SFA is failed")
        envirAssign("est.failed", TRUE)
        return(NULL)
    }
    
    beta <- coefs$beta
    if (!noSpatLag){
        coefs$rho = tail(beta, n=1)
        names(coefs$rho) = "Rho"
        coefs$beta = head(beta, -1)
    }else{
        coefs$rho = 0
    }
    p = olsen.reparam(coefs)
    result = c(p$gamma,p$nu,p$lambda,p$rho)
    logging("Initial values:", result, level="info")
    return(result)
}

spfrontier100HN.prepare = function(formula, data,W,...){
    prepareXY(formula, data)
    envirAssign("W", W)
}
spfrontier100HN.logL.direct <- function(parameters){
    k = length(parameters)
    beta = parameters[1:(k-3)]
    sigmaV = parameters[k-2]
    sigmaU = parameters[k-1]
    rho = parameters[k]
    par <- olsen.reparam(list(beta = beta, sigmaV=sigmaV,sigmaU=sigmaU, rho=rho))
    pargamma <- c(par$gamma,par$nu,par$lambda,par$rho)
    return(spfrontier100HN.logL(pargamma, quiet = T))
}

spfrontier100HN.postEstimation <- function(estimates){
    logging("spfrontier100HN: Handling estimates")
    if (status(estimates) == 0){
        coef <- olsen.reparamBack(spfrontier100HN.params(resultParams(estimates)))
        coefficients(estimates) <- coef
        hess <- optimHess(c(coef$beta,coef$sigmaV,coef$sigmaU,coef$rho),spfrontier100HN.logL.direct)
        hessian(estimates) <- hess
        
        X = envirGet("X")
        W = envirGet("W")
        y = envirGet("y")
        n <- length(y)
        spDet <- solve(diag(n) - coef$rho * W)
        fittedY <- spDet %*% X %*% coef$beta
        rownames(fittedY) <- rownames(X)
        colnames(fittedY) <- c("fittedted values")
        fitted(estimates) <- fittedY
        
        estimates <- frontierHN.handle.residuals(estimates)
        
    }
    logging("spfrontier100HN: Completed")
    return(estimates)
}