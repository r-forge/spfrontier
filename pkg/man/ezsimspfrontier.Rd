\name{ezsimspfrontier}
\alias{ezsimspfrontier}
\title{Spatial stochastic frontier model simulation tests}
\usage{
ezsimspfrontier(runs, autoSave = 0, params = params000, seed = NULL,
  inefficiency = "half-normal", logging = "info")
}
\arguments{
  \item{runs}{a number of simulated samples}

  \item{params}{a set with parameters to be used in
  simulation. There are predefined parameter sets:\cr
  params000 - a non-spatial stochastic frontier with
  half-normal inefficiencies\cr params000T - a non-spatial
  stochastic frontier with truncated normal
  inefficiencies\cr params100 - a stochastic frontier with
  spatial lags of a dependent variable and with half-normal
  inefficiencies\cr params100T - a stochastic frontier with
  spatial lags of a dependent variable and with truncated
  normal inefficiencies\cr params110 - stochastic frontier
  with spatial lags of a dependent variable and of a
  symmetric error component and with half-normal
  inefficiencies\cr params111 - stochastic frontier with
  spatial lags of a dependent variable, a symmetric error
  component, and an inefficiency error component and with
  half-normal inefficiencies\cr params011 - stochastic
  frontier with spatial lags of a symmetric error component
  and an inefficiency error component and with half-normal
  inefficiencies\cr params001 - stochastic frontier with
  spatial lags an inefficiency error component and with
  half-normal inefficiencies\cr}

  \item{autoSave}{save intermediate results to files. See
  \code{\link{ezsim}} for details.}

  \item{seed}{a state for random number generation in R. If
  NULL (default), the initial state is random. See
  \code{\link{set.seed}} for details.}

  \item{inefficiency}{sets the distribution for
  inefficiency error component. Possible values are
  'half-normal' (for half-normal distribution) and
  'truncated' (for truncated normal distribution). By
  default set to 'half-normal'. See references for
  explanations}

  \item{logging}{an optional level of logging. Possible
  values are 'quiet','warn','info','debug'. By default set
  to quiet.}
}
\description{
\code{ezsimspfrontier} tests estimators of a spatial
stochastic frontier model with different parameters
}
\details{
The \code{ezsimspfrontier} function executes multiple calls
of the \code{spfrontier} estimator on a simulated data set,
generated on the base of provided parameters. The resulting
estimates can be analysed for biasedness, efficiency, etc.
}
\examples{
#res <- ezsimspfrontier(100, params = params000,  seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(100, params = params000T, seed = 999, inefficiency = "truncated",logging = "info")
#res <- ezsimspfrontier(100, params = params100,  seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(100, params = params100T, seed = 999, inefficiency = "truncated",logging = "info")
#All tests above work as appropriate
#res <- ezsimspfrontier(10, params = params010, seed = 999, inefficiency = "half-normal",logging = "debug")
#A problem with sigmaV
#res <- ezsimspfrontier(10, params = params001, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params110, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params111, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params011, seed = 999, inefficiency = "half-normal",logging = "info")
}
\seealso{
\code{\link{ezsim}}
}
\keyword{frontier,}
\keyword{simulation}
\keyword{spatial}
\keyword{stochastic}

