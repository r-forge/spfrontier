

spfrontier.test <- function(){
    
    #ezsimspfrontier(10, params=params001, seed=0)
    #Сделать gridSearch
    
    data( airports2011 )
    
    W <- 1/as.matrix(dist(cbind(airports$lon, airports$lat)))
    colnames(W) <- airports$ICAO_code
    rownames(W) <- airports$ICAO_code
    W[which(W==Inf)] <- 0
    
    
    formula <- log(PAX) ~ log(runways) + log(checkins) +log (gates)
    #formula <- log(PAX) ~ log(runways) + log(checkins)
    
    ols <- lm(formula , data=airports)
    summary(ols )
    plot(density(residuals(ols)))
    skewness(residuals(ols))
    
    sfa <- sfa(formula , data=airports)
    summary(sfa )
    
    
    model <- spfrontier(formula , data=airports, logging="info",control=list(maxit=1000))
    summary(model )
    
    
    model <- spfrontier(formula , data=airports, W_y=W, logging="info",control=list(maxit=1000))
    summary(model )
    
    #11.17009670  0.17131972  0.90283885  0.19531389 -0.00135427  0.43227078  0.48969254 -0.08249509
    model <- spfrontier(formula , data=airports, W_y=W, W_v=W, logging="debug",control=list())
    summary(model )
    
    model <- spfrontier(formula , data=airports, W_v=W, logging="debug",control=list())
    summary(model )
    
    #Посчиталось - 1 нот дефинит
    model <- spfrontier(formula , data=airports, W_u=W, logging="debug",control=list())
    summary(model )
    
    
    
    
}

spfrontier.test2 <- function(){
    
    data( airports.spain.2010 )
    
    #formula <- log(revenue) ~ log(APM) + log(ATM) + log(staff_cost)
    #formula <- -log(ATM) ~ log(APM/ATM) + log(DA) + log(staff_cost)+ log(terminals)+ log(runways)
    
    formula <- -log(ATM) ~ log(APM/ATM) + log(DA) + log(staff_cost)
    
    data.spain2010 <- data.spain2010[-which.max(residuals(ols)),]
    
    ols <- lm(formula , data=data.spain2010)
    data.spain2010 <- data.spain2010[-which.max(residuals(ols)),]
    ols <- lm(formula , data=data.spain2010)
    summary(ols )
    plot(density(residuals(ols)))
    skewness(residuals(ols))
    
    
    sfa <- sfa(formula , data=data.spain2010)
    summary(sfa )
    plot(density(residuals(sfa )))
    
    W <- 1/as.matrix(dist(cbind(data.spain2010$lon, data.spain2010$lat)))
    colnames(W) <- data.spain2010$ICAO_code
    rownames(W) <- data.spain2010$ICAO_code
    W[which(W==Inf)] <- 0
    
    
    
    model <- spfrontier(formula , data=data.spain2010, logging="info",control=list(maxit=1000,reltol=1e-16))
    summary(model )
    
    model <- spfrontier(formula , data=data.spain2010, W_y=W, logging="info",control=list(maxit=1000,reltol=1e-16))
    summary(model )
    plot(density(residuals(model )))
    
    #Эта работает!formula <- -log(ATM) ~ log(APM/ATM) + log(DA) + log(staff_cost) - при удалении выброса (чтобы скос был отрицательным)
    model <- spfrontier(formula , data=data.spain2010, W_v=W, logging="info")
    summary(model )
    
    model <- spfrontier(formula , data=data.spain2010, W_y=W, W_v=W, logging="info")
    summary(model )
    
    model <- spfrontier(formula , data=data.spain2010,  W_u=W, logging="info",control=list())
    summary(model )
    
    
    model <- spfrontier(formula , data=data.spain2010, W_y=W, W_v=W,  W_u=W, logging="info")
    summary(model )
}


spfrontier.test3 <- function(){
    data(front41Data)
    cobbDouglas <- sfa( log( output ) ~ log( capital ) + log( labour ), data = front41Data )
    print(summary( cobbDouglas ))
    
    
    cobbDouglasT <- sfa( log( output ) ~ log( capital ) + log( labour ), data = front41Data, truncNorm=T )
    print(summary( cobbDouglasT ))
    
    
    model <- spfrontier( log( output ) ~ log( capital ) + log( labour ), data = front41Data, logging="info")
    print(summary( model ))
    
    model <- spfrontier( log( output ) ~ log( capital ) + log( labour ), data = front41Data, inefficiency="truncated",logging="info")
    print(summary( model ))
    
    
    W <- genW(nrow(front41Data))
    Wq <- genW(nrow(front41Data), type="queen")
    
    model <- spfrontier( log( output ) ~ log( capital ) + log( labour ), data = front41Data, W_y = W, logging="info")
    print(summary( model ))
    
    model <- spfrontier(log( output ) ~ log( capital ) + log( labour ), data = front41Data, W_y = W, W_v = Wq, logging="debug")
    print(summary(model))
    
}


spfrontier.test4 <- function(){
    
    data( airports.greece )
    formula <- log(WLU) ~ log(openning_hours) + log(runway_area) + log(terminal_area) +log(parking_area)
    
    sfa <- sfa(formula , data=airports.greece)
    summary(sfa )
    
    W <- 1/as.matrix(dist(cbind(airports.greece$lon, airports.greece$lat)))
    colnames(W) <- airports.greece$ICAO
    rownames(W) <- airports.greece$ICAO
    W[which(W==Inf)] <- 0
    
    
    
    
    model <- spfrontier(formula , data=airports.greece, logging="info")
    summary(model )
    
    
    model <- spfrontier(formula , data=airports.greece, logging="info", W_y=W)
    summary(model )
    
    # ini <- c(-5.624516608,1.894072073,-0.330861597,0.485801714,0.231013146,0.610346959,0.007689081,-0.027243936)
    model <- spfrontier(formula , data=airports.greece, logging="debug", W_v=W,initialValues=ini)
    summary(model )
    
    # -5.6234145  1.9277522 -0.3486990  0.4720372  0.2356252  0.5312435  1.0359464 -0.8000000
    model <- spfrontier(formula , data=airports.greece, logging="debug", W_u=W,control=list(grid.rhoU = 100))
    summary(model )
    
    # -4.997903703  1.887882894 -0.264050123  0.480395710  0.171142787 -0.005351067  0.461776359  1.373579604  0.000229417 -0.450620260 
    model <- spfrontier(formula , data=airports.greece, logging="debug", W_y=W, W_v=W, W_u=W)
    summary(model )
}