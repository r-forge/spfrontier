spfrontier111TN.params = function(parameters){
    k = length(parameters)
    p = spfrontier110HN.params(head(parameters,-2))
    p$rho3 = parameters[k-1]
    p$mu = parameters[k]
    return(p)
}

spfrontier111TN.logL <- function(parameters){
    counter = envirCounter("spfrontier111TN.logL")
    logging(paste("Evaluating 'spfrontier111TN.logL' for parameters","[run=",counter,"]:"),parameters) 
    y = envirGet("y")
    X = envirGet("X")
    W = envirGet("W")
    W2 = envirGet("W2")
    W3 = envirGet("W3")
    p = spfrontier111TN.params(parameters)
    ret = -10e8
    n = length(y)
    I <- diag(n)
    SpDet <- det(I-p$rho*W)
    
    if ((p$sigmaV>0) && (p$sigmaU>0) && (abs(p$rho)<1) && (abs(p$rho2)<1) && (abs(p$rho3)<1) && (SpDet>0)){
            
            e <- y    - p$rho * W %*% y -    X %*% p$beta
            
            Sp2 <- solve(I-p$rho2*W2)
            Sp3 <- solve(I-p$rho3*W3)
            
            mSigma = p$sigmaV^2*t(Sp2)%*%Sp2
            mOmega = p$sigmaU^2*t(Sp3)%*%Sp3
            mC = mSigma + mOmega
            mB = mOmega%*%solve(mC)%*%mSigma 
            mA = mB %*% solve(mSigma)
            mD = -mB %*% solve(mOmega)
            vMu = rep(p$mu, n)
            tryCatch({
                ret<- log(SpDet)-log(ptmvnorm(lowerx=rep(-Inf, n),upperx = rep(0, n),mean=as.vector(t(vMu)), sigma=mOmega))+log(dmvnorm(x=as.vector(e+vMu),mean=rep(0, n), sigma=mC))+log(ptmvnorm(lowerx=rep(-Inf, n),upperx = rep(0, n),mean=as.vector(t(-mA%*%e-mD%*%vMu)), sigma=mB))
            }, error = function(e){
                logging(e$message, level="warn")
            })
    }else{
        logging("Parameters are out of space")
    }
    if (is.nan(ret) || (ret==-Inf)) ret = -10e8
    logging(paste("spfrontier111TN.logL =",-ret,"\n"))
    return(-ret)
}

#sararsfa111.hnormal.lf <- function(parameters){
#    return(sararsfa111.tnormal.lf(c(parameters,0)))
#}
#sararsfa111.hnormal.ini<-function(formula, data){
#    ini = sararsfa111.tnormal.ini(formula, data)
#    return(head(ini, -1))
#}

spfrontier111TN.ini<-function(formula, data){
    logging("spfrontier111TN: calculating initial values")
    W = envirGet("W")
    iniModel <- spfrontier(formula,data,model="spfrontier100HN",W_y=W)
    coefs <- iniModel@coefficients
    if (status(iniModel) > 0){
        iniModel <- spfrontier(formula,data,model="frontierHN")
        coefs <- iniModel@coefficients
        coefs$rho <- 0
    }
    
    y = envirGet("y")
    X = envirGet("X")
    W2 = envirGet("W2")
    e = y - coefs$rho * y - X %*% coefs$beta
    W2e = W2 %*% e
    model = lm(e ~ W2e-1, data = data.frame(e, W2e))

    
    rho2 = coef(model)
    rho3 = 0
    ve = resid(model)
    mu = - mean(ve)
    result <-c(coefs$beta,coefs$sigmaV,coefs$sigmaU,coefs$rho, rho2, rho3, mu)
    names(result) <-c(paste("Beta", seq(from=0, to=length(coefs$beta)-1),sep = ""),"SigmaV","SigmaU", "Rho", "Rho2", "Rho3", "Mu")
    logging("Initial values:", result, level="info")
    return(result)
}

spfrontier111TN.prepare = function(formula, data,W,W2,W3,...){
    prepareXY(formula, data)
    envirAssign("W", W)
    envirAssign("W2", W2)
    envirAssign("W3", W3)
}

spfrontier111TN.postEstimation <- function(estimates){
    logging("spfrontier111TN: Handling estimates")
    if (status(estimates) == 0){
        p <- spfrontier111TN.params(resultParams(estimates))
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
        colnames(fittedY) <- c("Fitted values")
        fitted(estimates) <- fittedY
        
        resid <- y - fittedY
        rownames(resid) <- rownames(y)
        colnames(resid) <- c("Residuals")
        residuals(estimates) <- resid
        
        Sp2 <- solve(diag(n)-p$rho2*W2)
        Sp3 <- solve(diag(n)-p$rho3*W3)
        
        mSigma = p$sigmaV^2*t(Sp2)%*%Sp2
        mOmega = p$sigmaU^2*t(Sp3)%*%Sp3
        mC = mSigma + mOmega
        mB = mOmega%*%solve(mC)%*%mSigma 
        mA = mB %*% solve(mSigma)
        mD = -mB %*% solve(mOmega)
        vMu = rep(p$mu, n)
        
        u <- mtmvnorm(lower= rep(0, n),mean=as.vector(t(-mA%*%resid-mD%*%vMu)), sigma=mB, doComputeVariance=FALSE)$tmean
        eff <- cbind(exp(-u))
        rownames(eff) <- rownames(y)
        colnames(eff) <- c("Efficiency values")
        efficiencies(estimates) <- eff
        return(estimates)
    }
    logging("spfrontier111TN: Completed")
    return(estimates)
}