# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Definitions of S4 classes and methods
#

#
# Generic definitions
#
setGeneric("id",function(object){standardGeneric ("id")})
setGeneric("estimators",function(object){standardGeneric ("estimators")})
setGeneric("estimators<-",function(object,value){standardGeneric("estimators<-")})
setGeneric("execute",function(estimator, formula, data,initialValues,...){standardGeneric("execute")})
setGeneric("resultParams",function(object){standardGeneric ("resultParams")})
setGeneric("coefficients<-",function(object,value){standardGeneric("coefficients<-")})
setGeneric("residuals<-",function(object,value){standardGeneric("residuals<-")})
setGeneric("hessian<-",function(object,value){standardGeneric("hessian<-")})
setGeneric("stdErrors",function(object){standardGeneric("stdErrors")})
setGeneric("stdErrors<-",function(object,value){standardGeneric("stdErrors<-")})
setGeneric("fitted<-",function(object,value){standardGeneric("fitted<-")})
setGeneric("efficiencies",function(object){standardGeneric ("efficiencies")})
setGeneric("efficiencies<-",function(object,value){standardGeneric("efficiencies<-")})
setGeneric("status",function(object){standardGeneric ("status")})

#
# Class unions
#
setClassUnion("functionOrNULL", c("function", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))


#
# Class Model
# Stores information about model specifications and associated estimators.
#

setClass("Model", 
         representation(
             id="character", 
             name="character",
             estimators = "list"
             )
         )


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

# Accessors for the Model class
setMethod("id","Model", 
    function(object) {
        return(object@id)
    }
)

setMethod("estimators","Model", 
    function(object) {
        return(object@estimators)
    }
)

# Setters for the Model class
setReplaceMethod("estimators","Model", 
    function(object, value){
        object@estimators <- value
        return(object)
    }
)

#
# Class Estimator
# Stores information about MLE estimator. 
# MLE estimator includes three main functions (log-likelihood function,
# initial values provider, and gradient)
# and two service functions (initialise and handle postEstimation)
#

setClass("Estimator", 
         representation(
             id="character",
             initialize = "function",
             initialValues = "function",
             logL = "function",
             gradient = "functionOrNULL",
             postEstimation = "function")
         )
setMethod("show", "Estimator",
    function(object){
        cat("Estimator[ID=",object@id,"]\n", sep="")
    }
)

# Accessors for the Estimator class
setMethod("id","Estimator", 
    function(object) {
        return(object@id)
    }
)

# Executes the estimation procedure
setMethod("execute","Estimator", 
    function(estimator, formula, data, initialValues, ...){
        initEnvir(...)
        logging("Estimator started")
        estimator@initialize(formula, data,...)
        if (is.null(initialValues)){
            initialValues <- estimator@initialValues(formula, data)    
        }
        estimates <- optimEstimator(formula, data, 
                                     estimator@logL, 
                                     ini = initialValues, 
                                     gr=estimator@gradient)
        ret <- estimator@postEstimation(estimates)
        finalizeEnvir()
        return(ret)
    }
)


#
# Class ModelEstimates
# Stores information about estimated model
#
setClass("ModelEstimates", 
         representation(
             coefficients = "listOrNULL", 
             resultParams = "vector", 
             status = "numeric", 
             logL = "numeric", 
             logLcalls = "vector", 
             hessian = "matrix", 
             stdErrors = "vector", 
             residuals = "matrix", 
             fitted = "matrix", 
             efficiencies = "matrix")
)

setMethod("show", "ModelEstimates",
          function(object){
              cat("ModelEstimates\n", sep="")
              cat("Coefficients:\n")
              print(coefficients(object))
          }
)


# Accessors for the Estimator class

setMethod("coefficients", "ModelEstimates",
          function(object){
              return(object@coefficients)
          }
)

setMethod("resultParams", "ModelEstimates",
          function(object){
              return(unlist(object@resultParams))
          }
)

setMethod("fitted", "ModelEstimates",
          function(object){
              return(object@fitted)
          }
)

setMethod("efficiencies", "ModelEstimates",
          function(object){
              return(object@efficiencies)
          }
)


setMethod("residuals", "ModelEstimates",
          function(object){
              return(object@residuals)
          }
)

setMethod("stdErrors","ModelEstimates", 
          function(object) {
              return(object@stdErrors)
          }
)


setMethod("status","ModelEstimates", 
          function(object) {
              return(object@status)
        }
)

# Setters for the ModelEstimates class
setReplaceMethod("coefficients","ModelEstimates", 
                 function(object, value){
                     object@coefficients <- value;
                     return(object)
                }
)

setReplaceMethod("fitted","ModelEstimates", 
                 function(object, value){
                     object@fitted <- value;
                     return(object)
                }
)

setReplaceMethod("efficiencies","ModelEstimates", 
                 function(object, value){
                     object@efficiencies <- value;
                     return(object)
                }
)

setReplaceMethod("residuals","ModelEstimates", 
                 function(object, value){
                     object@residuals <- value;
                     return(object)
                }
)

setReplaceMethod("hessian","ModelEstimates", 
                 function(object, value){
                     object@hessian <- value;
                     OI <- solve(object@hessian)
                     object@stdErrors <- sqrt(diag(OI))
                     return(object)
                }
)

setReplaceMethod("stdErrors","ModelEstimates", 
                 function(object, value){
                    object@stdErrors <- value;
                    return(object)
                }
)

# Summary of the estimated model
setMethod("summary", "ModelEstimates",
          function(object){
              cat("ModelEstimates\n", sep="")
              cat("Status: ", object@status,"\n")
              cat("Estimates:\n")
              coef <- unlist(coefficients(object))
              zval <- coef/object@stdErrors
              pval <- 2*(1-pnorm(abs(zval)))
              digits <- max(5, .Options$digits - 3)
              e <- cbind(format(signif(coef,digits=digits)), 
                         format(signif(object@stdErrors,digits=digits)), 
                         format(signif(zval,digits=digits)),
                         format.pval(pval,digits=digits),
                         sapply(pval,pvalMark))
              rownames(e) <- names(coef)
              colnames(e) <- c("Estimate","Std. Error", "z value", "Pr(>|z|)","")
              print(e,quote = FALSE)
              cat("---\n")
              cat("Signif. codes:    0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")
              cat("Log-likelihood function value: ", object@logL,"\n")
              cat("Log-likelihood function/gradient calls: ", unlist(object@logLcalls),"\n")
          }
)
