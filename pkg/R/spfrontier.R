spfrontier <- function(formula, data,
                       W_y = NULL, W_v = NULL,W_u = NULL,
                       model="frontierHN",
                       method="mle",
                       initialValues=NULL,
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
    
    
    modelEstimates <- execute(estimator, formula,data,initialValues=initialValues,W=W_y,W2=W_v,W3=W_u,loggingLevel=logging,control=control)
}

spfrontier.test <- function(){
    data( front41Data )
    print("Estimator front41Data")
    print(summary(sfa( log( output ) ~ log( capital ) + log( labour ), data = front41Data)))
    #debug((estimators(.Models[["spfrontier111TN"]])[[1]])@logL)
    #(estimators(.Models[["spfrontier110HN"]])[[1]])@gradient
    
    print("Model frontierHN")
    model <- spfrontier(log( output ) ~ log( capital ) + log( labour ), data = front41Data, logging="info", control=list())
    print(model)
    
    
    print("Model spfrontier100HN")
    W_y <- genW(nrow(front41Data))
    model <- spfrontier(log( output ) ~ log( capital ) + log( labour ),model="spfrontier100HN", W_y = W_y, data = front41Data, logging="quiet")
    print(summary(model))
    
    
    print("Model spfrontier110HN")
    W_y <- genW(nrow(front41Data))
    W_v <- genW(nrow(front41Data), type="queen")
    model <- spfrontier(log( output ) ~ log( capital ) + log( labour ),model="spfrontier110HN", W_y = W_y, W_v = W_v, data = front41Data, logging="debug")
    print(summary(model))
    
    print("Model spfrontier110HN")
    W_y <- genW(nrow(front41Data))
    W_v <- genW(nrow(front41Data), type="queen")
    W_u <- genW(nrow(front41Data), type="queen")
    model <- spfrontier(log( output ) ~ log( capital ) + log( labour ),model="spfrontier111TN", W_y = W_y, W_v = W_v, W_u = W_u, data = front41Data, logging="debug")
    print(summary(model))
}