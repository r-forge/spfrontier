spfrontier <- function(formula, data, 
                       W_y = NULL, W_v = NULL,W_u = NULL,  
                       model="frontierHN",
                       method="mle",
                       ini.values=NULL,
                       logging = "quiet",
                       control=NULL){
  #Validation of model parameters
  
  model <- .Models[[model]] 
  if (is.null(model)){
    stop("Specified model is not supported")
  }
  
  estimator <- estimators(model)[[method]]
  if (is.null(estimator)){
    stop("Specified estimator is not supported")
  }
  
  if (!(logging %in% c("quiet", "info", "debug"))){
    warning("Specified logging level is not supported. Logging is set to 'quiet'")
    logging = "quiet"
  }
  #End of model parameters' validation 
  
  
  modelEstimates <- run(estimator, formula,data,ini.values=ini.values,W=W_y,W2=W_v,W3=W_u,logging.level=logging,control=control)
}

test <- function(){
  data( front41Data )
  model <- spfrontier(log( output ) ~ log( capital ) + log( labour ), data = front41Data, logging="debug")
  coef(model)
  
  debug((estimators(.Models[["frontierHN"]])[[1]])@initialize)
  sfa( log( output ) ~ log( capital ) + log( labour ), data = front41Data)
}