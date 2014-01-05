.Models <- list()
.onLoad <- function(libname, pkgname){
  registerModel(new("Model", id = "frontierHN", 
                    name = "Non-spatial stohastic frontier model with half-normal inefficiencies"))
  registerModel(new("Model", id = "frontierTN", 
                    name = "Non-spatial stohastic frontier model with truncated normal inefficiencies"))
  registerModel(new("Model", id = "spfrontier100HN", 
                    name = "Spatial autoregressive stohastic frontier model with non-spatial disturbances and non-spatial half-normal inefficiencies"))
  registerModel(new("Model", id = "spfrontier110HN", 
                    name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and non-spatial half-normal inefficiencies"))
  registerModel(new("Model", id = "spfrontier111HN", 
                    name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and spatial autoregressive half-normal inefficiencies"))
  registerModel(new("Model", id = "spfrontier111TN", 
                    name = "Spatial autoregressive stohastic frontier model with spatial autoregressive disturbances and spatial autoregressive truncated normal inefficiencies"))
  
  registerEstimator("frontierHN",
                    new("Estimator", 
                        id = "mle", 
                        initialize = frontierHN.prepare, 
                        ini.values = frontierHN.ini, 
                        logL = frontierHN.logL, 
                        gradient = frontierHN.logL.gradient, 
                        handle.estimates = frontierHN.handle.estimates
                    ))
  
  registerEstimator("spfrontier100HN",
                    new("Estimator", 
                        id = "mle", 
                        initialize = spfrontier100HN.prepare, 
                        ini.values = spfrontier100HN.ini, 
                        logL = spfrontier100HN.logL, 
                        gradient = spfrontier100HN.logL.gradient, 
                        handle.estimates = spfrontier100HN.handle.estimates
                    ))
  
  registerEstimator("spfrontier110HN",
                    new("Estimator", 
                        id = "mle", 
                        initialize = spfrontier110HN.prepare, 
                        ini.values = spfrontier110HN.ini, 
                        logL = spfrontier110HN.logL,
                        handle.estimates = spfrontier110HN.handle.estimates
                    ))
  
  registerEstimator("spfrontier111TN",
                    new("Estimator", 
                        id = "mle", 
                        initialize = spfrontier111TN.prepare, 
                        ini.values = spfrontier111TN.ini, 
                        logL = spfrontier111TN.logL,
                        handle.estimates = spfrontier111TN.handle.estimates
                    ))
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage("Thank you for using spfrontier!")
}