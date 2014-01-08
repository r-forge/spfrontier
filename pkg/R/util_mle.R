# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Routines for maximum likelihood estimator
#

#
# Wrapper for 'optim' procedure from 'stats' package
#
optimEstimator <- function(formula, data, func.lf, ini,gr=NULL){
    ret <- prepareEstimates(status = 1000)
    if (is.null(ini)){
        logging("Ini is not defined",caller = match.call())
        return(prepareEstimates(status = 1001))
    }
    tryCatch({
        control <- envirGet("control")
        if(!is.null(control)) {
            logging("Optimizer control",control)
        }
        p<-optim(par=ini,fn=func.lf,gr=gr,
                 method="BFGS", control=control,hessian = T)
        logging(paste("Convergence is",
                      ifelse(p$convergence==0, "", "NOT"),"archieved",
                      "[",p$convergence,"]"),
                level = "info",
                caller = match.call())    
        logging(paste("Log-Likelihood value = ",p$value),
                level = "info",
                caller = match.call())  
        ret <- prepareEstimates(estimates = p)
    }, error = function(e){ 
        logging(e$message, level="warn")
    })
    return(ret)
}

#
# Creates the ModelEstimates class on the base of raw (optim) estimation results
#
prepareEstimates <- function(estimates = NULL, status = 0){
    ret <- new("ModelEstimates", status = status)
    if (!is.null(estimates)) {
        ret <- new("ModelEstimates",
                   resultParams = estimates$par,
                   status = estimates$convergence,
                   logL = estimates$value,
                   logLcalls = estimates$counts
                   )
        
        #Setting hessian and implicitly calculating standard errors
        if (!is.null(estimates$hessian)){
            hessian(ret) <- estimates$hessian
        }
    }
    return(ret)
}