.Models = list()

setClass("Model", representation(id="character", name="character" ,estimators = "list"))
setMethod("show", "Model",
  function(object){
    cat("Model[ID=",object@id,", name=",object@name,"]\n", sep="")
    if (length(object@estimators) == 0) {
      cat("No estimators are registered\n", sep="")
    }else{
      cat("Registered estimators:\n", sep="")
      for (e in object@estimators) {
        show(e)
      }
    }
  }
)

setGenericVerif <- function(x,y){if(!isGeneric(x)){setGeneric(x,y)}else{}}

setGenericVerif("getId",function(object){standardGeneric ("getId")})
setMethod("getId","Model", function(object){return(object@id)})

setGenericVerif("estimators",function(model){standardGeneric ("estimators")})
setGenericVerif("estimators<-",function(model,value){standardGeneric("estimators<-")})

setMethod("estimators","Model", function(model){return(model@estimators)})
setReplaceMethod("estimators","Model", function(model, value){
  model@estimators <- value;
  return(model)
})


setClass("Estimator", representation(id="character" ,
                                     initialize = "function", 
                                     ini.values = "function", 
                                     logL = "function", 
                                     gradient = "function", 
                                     handle.estimates = "function"))
setMethod("show", "Estimator",
    function(object){
      cat("Estimator[ID=",object@id,"]\n", sep="")
    }
)
setMethod("getId","Estimator", function(object){return(object@id)})

setGenericVerif("run",function(estimator, formula, data,ini.values,...){standardGeneric("run")})
setMethod("run","Estimator", function(estimator, formula, data,ini.values, ...){
  envir.init(...)
  logger.start()
  logger.debug("Estimator started")
  estimator@initialize(formula, data,...)
  if (is.null(ini.values)){
    ini.values <- estimator@ini.values(formula, data)  
  } 
  estimates <- optim.estimator(formula, data, estimator@logL, ini = ini.values, gr=estimator@gradient)
  ret <- estimator@handle.estimates(estimates)
  envir.finalise()
  return(ret)
})



setClass("ModelEstimates", representation(
                                     coefficients = "list", 
                                     status = "numeric", 
                                     logL = "numeric", 
                                     hessian = "matrix", 
                                     residuals = "list", 
                                     fit = "list"))
setMethod("show", "ModelEstimates",
          function(object){
            cat("ModelEstimates\n", sep="")
            print(object@coefficients)
          }
)
setMethod("coef", "ModelEstimates",
          function(object){
            return(unlist(object@coefficients))
          }
)
setGenericVerif("status",function(object){standardGeneric ("status")})
setMethod("status","ModelEstimates", function(object){return(object@status)})

registerModel = function(model){
  .Models[[getId(model)]] <<- model
}


registerEstimator = function(modelId, estimator){
  model <- .Models[[modelId]]
  if (is.null(model)){
    stop(cat("No model is registered for ID=", modelId))
  }
  est <- estimators(model)
  if (!is.null(est[[getId(estimator)]]))
    warning(paste("Estimator",getId(estimator)," is already registered for the model",modelId," - replacing"))
  est[[getId(estimator)]] <- estimator
  estimators(model) <- est
  registerModel(model)
}


registerModel(new("Model", id = "frontierHN", name = "Non-spatial stohastic frontier model with half-normal inefficiencies"))
registerModel(new("Model", id = "frontierTN", name = "Non-spatial stohastic frontier model with truncated normal inefficiencies"))
registerModel(new("Model", id = "spfrontier100HN", name = "Spatial autoregressive stohastic frontier model with non-spatial disturbances and non-spatial half-normal inefficiencies"))
registerModel(new("Model", id = "spfrontier101HN", name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and non-spatial half-normal inefficiencies"))
registerModel(new("Model", id = "spfrontier111HN", name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and spatial autoregressive half-normal inefficiencies"))
registerModel(new("Model", id = "spfrontier111TN", name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and spatial autoregressive truncated normal inefficiencies"))

.Models