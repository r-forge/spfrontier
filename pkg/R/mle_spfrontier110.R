spfrontier110HN.params = function(parameters){
    k = length(parameters)
    beta <- parameters[1:(k-4)]
    sigmaV <- parameters[k-3]
    sigmaU <- parameters[k-2]
    rho <- parameters[k-1]
    rho2 <- parameters[k]
    
    
    X = envirGet("X")
    names(beta) = colnames(X)
    names(sigmaV) = "SigmaV"
    names(sigmaU) = "SigmaU"
    names(rho) = "Rho"
    names(rho2) = "Rho2"
    
    return(list(beta=beta,sigmaV=sigmaV,sigmaU=sigmaU,rho=rho,rho2=rho2))
}

spfrontier110HN.logL <- function(parameters){ 
    counter = envirCounter("spfrontier110HN.logL")
    logging(paste("Evaluating 'spfrontier110HN.logL' for parameters","[run=",counter,"]:"),parameters) 
    y = envirGet("y")
    X = envirGet("X")
    W = envirGet("W")
    W2 = envirGet("W2")
    p = spfrontier110HN.params(parameters)
    ret = -10e8
    if ((p$sigmaV>0) && (p$sigmaU>0) && (abs(p$rho)<1) && (abs(p$rho2)<1)){ 
        n = length(y)
        e <- y    - p$rho * W %*% y -    X %*% p$beta
        I <- diag(n)        
        Sp2 <- solve(I-p$rho2*W2)
        mCapitalSigma <- t(Sp2)%*%Sp2
        mCapitalTheta <- p$sigmaU^2*I+p$sigmaV^2*mCapitalSigma
        mCapitalOmega <-p$sigmaU^2*p$sigmaV^2*mCapitalSigma %*% solve(mCapitalTheta)
        mMu <- -p$sigmaV^(-2)*mCapitalOmega%*%solve(mCapitalSigma)%*%e
        ret<- log(det(I-p$rho*W))+log(ptmvnorm(lowerx=rep(0, n),upperx = rep( Inf, n),mean=as.vector(t(mMu)), sigma=mCapitalOmega))+log(dmvnorm(x=as.vector(e),mean=rep(0, n), sigma=mCapitalTheta))
    }else{
     logging("Parameters are out of space")
    }
    if (is.nan(ret) || (ret==-Inf)) ret = -10e8
    logging(paste("spfrontier110HN.logL =",-ret,"\n"))
    return(-ret)
}


spfrontier110HN.ini<-function(formula, data){
    logging("spfrontier110HN: calculating initial values")
    W = envirGet("W")
    iniModel <- spfrontier(formula,data,model="spfrontier100HN",W_y=W)
    coefs <- iniModel@coefficients
    if (status(iniModel) > 0){
        iniModel <- spfrontier(formula,data,model="frontierHN")
        coefs <- iniModel@coefficients
        coefs$rho <- 0
    }
    if (length(coefs$beta)==0) {
        print("coefs is failed")
        return(NULL)
    }
    rho2 = 0
    result <-c(coefs$beta,coefs$sigmaV,coefs$sigmaU,coefs$rho, rho2)
    names(result) <-c(paste("Beta", seq(from=0, to=length(coefs$beta)-1),sep = ""),"SigmaV","SigmaU", "Rho", "Rho2")
    logging("Initial values:", result, level="info")
    return(result)
}

spfrontier110HN.prepare = function(formula, data,W,W2,...){
    prepareXY(formula, data)
    envirAssign("W", W)
    envirAssign("W2", W2)
}

spfrontier110HN.postEstimation <- function(estimates){
    logging("spfrontier110HN: Handling estimates")
    if (status(estimates) == 0){
        p <- spfrontier110HN.params(resultParams(estimates))
        
        coefficients(estimates) <- p
        
        y = envirGet("y")
        X = envirGet("X")
        W = envirGet("W")
        W2 = envirGet("W2")
        W3 = envirGet("W3")
        
        n <- length(y)
        spDet <- solve(diag(n) - p$rho * W)
        fittedY <- spDet %*% X %*% p$beta
        rownames(fittedY) <- rownames(X)
        colnames(fittedY) <- c("fittedted values")
        fitted(estimates) <- fittedY
        
        resid <- y - fittedY
        rownames(resid) <- rownames(y)
        colnames(resid) <- c("Residuals")
        residuals(estimates) <- resid
        
        Sp2 <- solve(diag(n)-p$rho2*W2)
        Sp3 <- diag(n)
        
        mSigma = p$sigmaV^2*t(Sp2)%*%Sp2
        mOmega = p$sigmaU^2*t(Sp3)%*%Sp3
        mC = mSigma + mOmega
        mB = mOmega%*%solve(mC)%*%mSigma 
        mA = mB %*% solve(mSigma)
        mD = -mB %*% solve(mOmega)
        vMu = rep(0, n)
        
        u <- mtmvnorm(lower= rep(0, n),mean=as.vector(t(-mA%*%resid-mD%*%vMu)), sigma=mB, doComputeVariance=FALSE)$tmean
        eff <- cbind(exp(-u))
        rownames(eff) <- rownames(y)
        colnames(eff) <- c("Efficiency values")
        efficiencies(estimates) <- eff
        return(estimates)
    }
    logging("spfrontier110HN: Completed")
    return(estimates)
}