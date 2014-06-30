params <- list(n=100,
                  beta0=5,
                  beta1=10,
                  beta2=1,
                  sigmaV=0.5, 
                  sigmaU=2.5,
               rhoU=0.2)

parDef <- createParDef(selection = params, banker = list(loggingLevel="debug",inefficiency="half-normal",control=list(),parDef=createParDef(selection = params, banker=list())))
set.seed(999)
dgp <- evalFunctionOnParameterDef(parDef, spfrontier.dgp)

values <- dgp$tv
logLikelihood(formula=dgp$formula, data=dgp$data,
              W_y = dgp$W_y, W_v = dgp$W_v,W_u = dgp$W_u,
              inefficiency = dgp$inefficiency,
              values=values,
              logging = "debug")