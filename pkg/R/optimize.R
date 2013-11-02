
optim.estimator <- function(formula, data, func.lf, func.ini,gr=NULL,silent=TRUE, ...){
  ini <- func.ini(formula, data)
  return(optim.estimator.withIni(formula, data, func.lf, ini,gr,...,silent))
}

optim.estimator.withIni <- function(formula, data, func.lf, ini,gr=NULL,silent=TRUE, ...){
  result <- NULL
  mf <- model.frame(formula, data)
  
  
  y <- as.matrix(model.response(mf))
  X <- as.matrix(mf[-1])
 
  tm <- attr(mf, "terms")
  intercept <- attr(tm, "intercept") == 1
  if (intercept)  X <- cbind(1L,X)
  tryCatch({
    p<-optim(ini,func.lf,method="BFGS",hessian=T, y, X, gr=gr, control=list(maxit=10000,reltol=1e-12))
    if (!silent){
      print("Estimates: ")
      print(p$par)
      print(paste("Convergence: ",p$convergence))
      print(paste("LogLikelihood value: ",p$value))
      print("----------------------------------------")
      flush.console()
    }
    if (p$convergence==0) result <-list(estimate=p$par,value=p$value)
  }, error = function(e) if (!silent){print(e$message)})
  return(result)
}