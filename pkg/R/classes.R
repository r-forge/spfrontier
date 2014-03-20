# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Definitions of S4 classes and methods
#

#
# Generic definitions
#
setGeneric("resultParams",function(object){standardGeneric ("resultParams")})
setGeneric("coefficients<-",function(object,value){standardGeneric("coefficients<-")})
setGeneric("residuals<-",function(object,value){standardGeneric("residuals<-")})
setGeneric("hessian",function(object){standardGeneric("hessian")})
setGeneric("hessian<-",function(object,value){standardGeneric("hessian<-")})
setGeneric("stdErrors",function(object){standardGeneric("stdErrors")})
setGeneric("stdErrors<-",function(object,value){standardGeneric("stdErrors<-")})
setGeneric("fitted<-",function(object,value){standardGeneric("fitted<-")})
setGeneric("efficiencies",function(object){standardGeneric ("efficiencies")})
setGeneric("efficiencies<-",function(object,value){standardGeneric("efficiencies<-")})
setGeneric("status",function(object){standardGeneric ("status")})

#
# Class ModelEstimates
# Stores information about estimated model
#
setClass("ModelEstimates", 
         representation(
             coefficients = "vector", 
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
setMethod("hessian","ModelEstimates", 
           function(object) {
               return(object@hessian)
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
                     tryCatch({
                         OI <- solve(object@hessian)
                         v <- diag(OI)
                         se <- c()
                         for (i in v){
                             if (i>0){
                                 se <- c(se, sqrt(i))
                             } else{
                                 se <- c(se, NA)
                             }
                         }
                         object@stdErrors <- se
                     }, error = function(e){
                         logging(e$message, level="warn")
                     })
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
              if (object@status == 0){
                  cat("Estimates:\n")
                  coef <- coefficients(object)
                  coef <- c(coef$beta,coef$rhoY, coef$sigmaV, coef$sigmaU, coef$rhoV, coef$rhoU, coef$mu)
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
          }
)
