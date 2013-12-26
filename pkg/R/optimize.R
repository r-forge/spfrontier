
optim.estimator <- function(formula, data, func.lf, func.ini,env,gr=NULL,silent=TRUE){
  ini = NULL
  if (exists("ini_value", envir =env)){
    ini = get("ini_value", envir =env)
  }else{
    ini <- func.ini(formula, data, env=env)  
  } 
  return(optim.estimator.withIni(formula, data, func.lf, ini,env=env,gr=gr,silent=silent))
}

optim.estimator.withIni <- function(formula, data, func.lf, ini,env=env,gr=NULL,silent=TRUE){
  result <- NULL
  mf <- model.frame(formula, data)
  
  
  y <- as.matrix(model.response(mf))
  X <- as.matrix(mf[-1])
  tm <- attr(mf, "terms")
  intercept <- attr(tm, "intercept") == 1
  if (intercept)  X <- cbind(1L,X)
  assign("X", X, envir=env)
  assign("y", y, envir=env)
  tryCatch({
    if (!silent){
      print("----------------------------------------")
      print("INI: ")
      print(ini)
    }
    if (is.null(ini)){
      stop("Ini is not defined")
    }
    p<-optim(par=ini,fn=func.lf,env=env,gr=gr,method="BFGS", control=list(maxit=10000,reltol=1e-16),hessian=T)
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