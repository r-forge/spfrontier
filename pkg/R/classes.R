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

setGenericVerif("estimators",function(object){standardGeneric ("estimators")})
setGenericVerif("estimators<-",function(object,value){standardGeneric("estimators<-")})

setMethod("estimators","Model", function(object){return(object@estimators)})
setReplaceMethod("estimators","Model", function(object, value){
  object@estimators <- value;
  return(object)
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
setGenericVerif("run",function(object, formula, data,ini.values,...){standardGeneric("run")})
setMethod("run","Estimator", function(object, formula, data,ini.values, ...){
  envir.init(...)
  logger.start()
  logger.debug("Estimator started")
  print(object@initialize)
  object@initialize(formula, data,...)
  if (is.null(ini.values)){
    ini.values <- object@ini.values(formula, data)  
  } 
  estimates <- optim.estimator(formula, data, object@logL, ini = ini.values, gr=object@gradient)
  ret <- object@handle.estimates(estimates)
  envir.finalise()
  return(ret)
})



setClass("ModelEstimates", representation(
                                     coefficients = "list", 
                                     status = "numeric", 
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

registerModel = function(model){
  .Models[[getId(model)]] <<- model
}

#Fix @
registerEstimator = function(modelId, estimator){
  model <- .Models[[modelId]]
  if (is.null(model)){
    stop(cat("No model is registered for ID=", modelId))
  }
  est <- estimators(model)
  if (!is.null(est[[estimator@id]]))
    warning(paste("Estimator",estimator@id," is already registered for the model",modelId," - replacing"))
  est[[estimator@id]] <- estimator
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