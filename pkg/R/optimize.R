optim.estimator <- function(formula, data, func.lf, ini,gr=NULL){
  result <- NULL
  if (is.null(ini)){
    logger.debug("Ini is not defined",caller = match.call())
    stop("Ini is not defined")
  }
  tryCatch({
    p<-optim(par=ini,fn=func.lf,gr=gr,method="BFGS", control=list(maxit=10000,reltol=1e-16),hessian=T)
    logger.info(paste("Convergence is",ifelse(p$convergence==0, "", "NOT"),"archieved","(",p$convergence,").","Log Likelihood value: ",p$value)
                ,caller = match.call())  
    if (p$convergence==0){
      result <-list(estimate=p$par,value=p$value)
    } 
  }, error = function(e) logger.warn(e$message))
  return(result)
}