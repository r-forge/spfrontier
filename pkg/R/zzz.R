# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Package initialization
#

# Global repository of available models
.Models <- list()

#
# Registers a model in the model repository 
#
registerModel = function(model){
    .Models[[id(model)]] <<- model
}

#
# Registers an estimator for a specified model
#
registerEstimator = function(modelId, estimator){
    model <- .Models[[modelId]]
    if (is.null(model)){
        stop(cat("No model is registered for ID=", modelId))
    }
    est <- estimators(model)
    if (!is.null(est[[id(estimator)]]))
        warning(paste("Estimator",id(estimator)," is already registered for the model",modelId," - replacing"))
    est[[id(estimator)]] <- estimator
    estimators(model) <- est
    registerModel(model)
}

#
# Package loading
#
.onLoad <- function(libname, pkgname){
    
    # Register supported models
    
    registerModel(new("Model", 
                      id = "frontierHN",
                      name = "Non-spatial stohastic frontier model with half-normal inefficiencies"))
    registerModel(new("Model", 
                      id = "frontierTN",
                      name = "Non-spatial stohastic frontier model with truncated normal inefficiencies"))
    registerModel(new("Model", 
                      id = "spfrontier100HN",
                      name = "Spatial autoregressive stohastic frontier model with non-spatial disturbances and non-spatial half-normal inefficiencies"))
    registerModel(new("Model", 
                      id = "spfrontier110HN",
                      name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and non-spatial half-normal inefficiencies"))
    registerModel(new("Model", 
                      id = "spfrontier111HN",
                      name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and spatial autoregressive half-normal inefficiencies"))
    registerModel(new("Model", 
                      id = "spfrontier111TN",
                      name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and spatial autoregressive truncated normal inefficiencies"))
    
    # Register implemented estimators
    registerEstimator("frontierHN",
                            new("Estimator", 
                                    id = "mle", 
                                    initialize = frontierHN.prepare, 
                                    initialValues = frontierHN.ini, 
                                    logL = frontierHN.logL, 
                                    gradient = frontierHN.logL.gradient, 
                                    postEstimation = frontierHN.postEstimation
                            ))
    
    registerEstimator("spfrontier100HN",
                            new("Estimator", 
                                    id = "mle", 
                                    initialize = spfrontier100HN.prepare, 
                                    initialValues = spfrontier100HN.ini, 
                                    logL = spfrontier100HN.logL, 
                                    gradient = spfrontier100HN.logL.gradient, 
                                    postEstimation = spfrontier100HN.postEstimation
                            ))
    
    registerEstimator("spfrontier110HN",
                            new("Estimator", 
                                    id = "mle", 
                                    initialize = spfrontier110HN.prepare, 
                                    initialValues = spfrontier110HN.ini, 
                                    logL = spfrontier110HN.logL,
                                    postEstimation = spfrontier110HN.postEstimation
                            ))
    
    registerEstimator("spfrontier111TN",
                            new("Estimator", 
                                    id = "mle", 
                                    initialize = spfrontier111TN.prepare, 
                                    initialValues = spfrontier111TN.ini, 
                                    logL = spfrontier111TN.logL,
                                    postEstimation = spfrontier111TN.postEstimation
                            ))
}

#
# Package attaching
#
.onAttach <- function(libname, pkgname){
    packageStartupMessage("Thank you for using 'spfrontier'")
}