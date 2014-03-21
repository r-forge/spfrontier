\name{spfrontier}
\alias{spfrontier}
\title{Spatial stochastic frontier model}
\usage{
spfrontier(formula, data, W_y = NULL, W_v = NULL, W_u = NULL,
  inefficiency = "half-normal", initialValues = NULL, logging = c("quiet",
  "info", "debug"), control = NULL, onlyCoef = F)
}
\arguments{
  \item{formula}{an object of class
  "\code{\link{formula}}": a symbolic description of the
  model to be fitted. The details of model specification
  are given under 'Details'.}

  \item{data}{data frame, containing the variables in the
  model}

  \item{W_y}{a spatial weight matrix for spatial lag of the
  dependent variable}

  \item{W_v}{a spatial weight matrix for spatial lag of the
  symmetric error term}

  \item{W_u}{a spatial weight matrix for spatial lag of the
  inefficiency error term}

  \item{initialValues}{an optional vector of initial
  values, used by maximum likelihood estimator. If not
  defined, estimator-specific method of initial values
  estimation is used.}

  \item{logging}{an optional level of logging. Possible
  values are 'quiet','warn','info','debug'. By default set
  to quiet.}

  \item{inefficiency}{sets the distribution for
  inefficiency error component. Possible values are
  'half-normal' (for half-normal distribution) and
  'truncated' (for truncated normal distribution). By
  default set to 'half-normal'. See references for
  explanations}

  \item{onlyCoef}{allows calculating only estimates for
  coefficients (with inefficiencies and other additional
  statistics). Developed generally for testing, to speed up
  the process.}

  \item{control}{an optional list of control parameters,
  passed to \code{\link{optim}} estimator from the
  '\code{\link{stats}} package}
}
\description{
\code{spfrontier} estimates spatial specifications of the
stochastic frontier model.
}
\details{
Models for estimation are specified symbolically, but
without any spatial components. Spatial components are
included implicitly on the base of the \code{model}
argument.
}
\examples{
data( airports )
W <- constructW(cbind(airports$lon, airports$lat),airports$ICAO_code)

formula <- log(PAX) ~ log(runways) + log(checkins) +log (gates)
ols <- lm(formula , data=airports)
summary(ols )
plot(density(stats::residuals(ols)))
skewness(stats::residuals(ols))

model <- spfrontier(formula , data=airports)
summary(model )

model <- spfrontier(formula , data=airports, W_y=W)
summary(model )

# Takes a long time
#model <- spfrontier(formula , data=airports, W_y=W, W_v=W, logging="debug",control=list())
#summary(model )



data(airports.greece)
formula <- log(WLU) ~ log(openning_hours) + log(runway_area) + log(terminal_area)
W <- constructW(cbind(airports.greece$lon, airports.greece$lat),airports.greece$ICAO)



model000 <- spfrontier(formula , data=airports.greece, logging="info")
summary(model000 )

model100 <- spfrontier(formula , data=airports.greece, logging="info", W_y=W)
summary(model100 )

# Takes a long time
#model010 <- spfrontier(formula , data=airports.greece, logging="debug", W_v=W)
#summary(model010 )
#efficiencies(model010)

# Takes a long time
#model110 <- spfrontier(formula , data=airports.greece, logging="debug", W_y=W, W_v=W)
#summary(model110)
#efficiencies(model110)

# Takes a long time
#model001 <- spfrontier(formula , data=airports.greece, logging="debug", W_u=W)
#summary(model001)
#efficiencies(model001)



data(airports.spain)
formula <- -log(ATM) ~ log(APM/ATM) + log(DA) + log(staff_cost)
W <- constructW(cbind(airports.spain$lon, airports.spain$lat),airports.spain$ICAO_code)
}
\references{
Kumbhakar, S.C. and Lovell, C.A.K (2000), Stochastic
Frontier Analysis, Cambridge University Press, U.K.
}
\keyword{frontier}
\keyword{spatial}
\keyword{stochastic}

