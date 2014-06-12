params000 <- list(n=c(100,200,300),
                  sigmaX=10, 
                  beta0=1,
                  beta1=-2,
                  beta2=3, 
                  sigmaV=0.2, 
                  sigmaU=1)

params000T <- params000
params000T$mu <- 0.4

params100 <- params000
params100$rhoY <- 0.2

params100T <- params000T
params100T$rhoY <- 0.2


params110 <- params100
params110$rhoV <- 0.3

params101 <- params100
params110$rhoU <- 0.5

params010 <- params110
params010$rhoY <- NULL

params111 <- params110
params111$rhoU <- 0.5

params011 <- params111
params011$rhoY <- NULL

params001 <- params011
params001$rhoV <- NULL

res000 <- ezsimspfrontier(100, params = params000,  seed = 999, inefficiency = "half-normal",logging = "info")
res000T <- ezsimspfrontier(100, params = params000T, seed = 999, inefficiency = "truncated",logging = "info")

res100 <- ezsimspfrontier(100, params = params100,  seed = 999, inefficiency = "half-normal",logging = "info")
res100_bias <- ezsimspfrontier(100, params = params100,  seed = 999, inefficiency = "half-normal",logging = "info", control=list(ignoreWy=T))

res100T <- ezsimspfrontier(100, params = params100T, seed = 999, inefficiency = "truncated",logging = "info")

params001A <- params001
params001A$rhoY <- 0
res001A <- ezsimspfrontier(100, params = params001A,  seed = 999, inefficiency = "half-normal",logging = "info", control=list(ignoreWu=T,replaceWyWu=T))

params010A <- params010
params010A$rhoY <- 0
res010A <- ezsimspfrontier(100, params = params010A,  seed = 999, inefficiency = "half-normal",logging = "info", control=list(ignoreWv=T,replaceWyWv=T))




#All tests above work as appropriate

#res <- ezsimspfrontier(10, params = params010, seed = 999, inefficiency = "half-normal",logging = "debug")
#A problem with sigmaV
#res001 <- ezsimspfrontier(100, params = params001, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params110, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params101, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params011, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params111, seed = 999, inefficiency = "half-normal",logging = "info")

res <- res100_bias

summary(res)
#subset=list(estimator=c("Intercept","X1","rhoY","sigmaV","sigmaU","mu")),
ps <- plot(res, parameters_priority = c("beta0","beta1","beta2","sigmaV","sigmaU","mu","rhoY","rhoV","rhoU"), return_print =T)
do.call(grid.arrange,  ps)
ps <- plot(res, type="density",parameters_priority = c("Intercept","X1","X2","rhoY","sigmaV","sigmaU","mu","rhoV","rhoU"), return_print =T)
do.call(grid.arrange,  ps)

