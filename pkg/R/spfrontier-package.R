#' Spatial Stochastic Frontier 
#' 
#' A set of tools for estimation (MLE) of various spatial specifications of stochastic frontier models
#' @name spfrontier-package
#' @docType package
#' @title Spatial Stochastic Frontier 
#' @author Dmitry Pavlyuk \email{Dmitry.V.Pavlyuk@@gmail.com}
#' @keywords spatial stochastic frontier
NULL

#' European airports statistical data
#' 
#' The \code{spfrontier} package includes the dataset \code{airports}, 
#' containing information about European airports infrastructure and traffic statistics in 2011.
#' 
#' 
#' @name airports
#' @rdname data-airports
#' @docType data
#' @format An unbalanced panel of 395 Euripean airports in 2008-2012 (1763  observations) on the following 31 variables.
#'
#' \describe{
#' 
#' \item{ICAO}{Airport ICAO code}
#' \item{AirportName}{Airport official name}
#' \item{Country}{Airport's country name}
#' \item{longitude}{Airport longitude}
#' \item{latitude}{Airport latitude}
#' \item{Year}{Observation year}
#' \item{PAX}{A number of carried passengers}
#' \item{ATM}{A number of of air transport movements served by an airport}
#' \item{Cargo}{A total volume of cargo served by an airport}
#' \item{Population100km}{A number of inhabitants, living in 100 km around an airport}
#' \item{Population200km}{A number of inhabitants, living in 200 km around an airport}
#' \item{Island}{1 if an airport is located on an island; 0 otherwise}
#' \item{GDPpc}{Gross domestic product per capita in airport's NUTS3 region}
#' \item{RevenueTotal}{Airport total revenue}
#' \item{RevenueAviation}{Airport aviation revenue}
#' \item{RevenueNonAviation}{Airport non-aviation revenue}
#' \item{RevenueHandling}{Airport revenue from handling services}
#' \item{RevenueParking}{Airport revenue from parking services}
#' \item{EBITDA}{Airport earnings before interest, taxes, depreciation, and amortization}
#' \item{NetProfit}{Airport net profit}
#' \item{DA}{Airport depreciation, and amortization}
#' \item{StaffCount}{A number of staff employed by an airport}
#' \item{StaffCost}{Airport staff cost}
#' \item{RunwayCount}{A number of airport runways}
#' \item{CheckinCount}{A number of airport check-iun facilities}
#' \item{GateCount}{A number of airport gates}
#' \item{TerminalCount}{A number of airport terminals}
#' \item{ParkingSpaces}{A number of airport parking spaces}
#' \item{RoutesDeparture}{A number of departure routes, served by an airport}
#' \item{RoutesArrival}{A number of arrival routes, served by an airport}
#' \item{Routes}{(RoutesDeparture + RoutesArrival)/2}
#' 
#' }
#' @source 
#' \describe{
#' \item{eurostat_report}{Eurostat (2013). European Statistics Database, Statistical Office of the European Communities (Eurostat)}
#' \item{report}{Airports' statistical reports(2011)}
#' \item{openflights}{Open Flights: Airport, airline and route data http://openflights.org/ (2013-05-31)}
#' \item{aena}{TDC (2012). Informe de fiscalizacion de la imputacion por la entidad "Aeropuertos Espanoles y Navegacion Aerea" (AENA) a cada uno de los aeropuertos de los ingresos, gastos, e inversiones correspondientes a la actividad aeroportuaria, en los ejercicios 2009 y 2010., Tribunal de Cuentas, Spain, Doc 938.}
#' \item{ciesin}{CIESIN, Columbia University. Gridded Population of the World: Future Estimates (GPWFE). (2005)}
#' }
NULL



#' Greece airports statistical data
#' 
#' The \code{spfrontier} package includes the dataset \code{airports}, 
#' containing information about Greece airports infrastructure and traffic statistics in 2011.
#' 
#' 
#' @name airports.greece
#' @rdname data-airports-greece
#' @docType data
#' @format A dataframe with 39 observations on the following 24 variables.
#'
#' \describe{
#' 
#' \item{name}{Airport title}
#' \item{ICAO_code}{Airport ICAO code}
#' \item{lat}{Airport latitude}
#' \item{lon}{Airport longitude}
#' \item{APM_winter}{A number of passengers carried during winter period}
#' \item{APM_summer}{A number of passengers carried during summer period}
#' \item{APM}{A number of passengers carried (winter + summer)}
#' \item{cargo_winter}{A total volume of cargo served by an airport during winter period}
#' \item{cargo_summer}{A total volume of cargo served by an airport during summer period}
#' \item{cargo}{A number volume of cargo served by an airport (winter + summer)}
#' \item{ATM_winter}{A number of air transport movements served by an airport during winter period}
#' \item{ATM_summer}{A number of air transport movements served by an airport during summer period}
#' \item{ATM}{A number of air transport movements served by an airport (winter + summer)}
#' \item{openning_hours_winter}{A total number openning hours during winter period}
#' \item{openning_hours_summer}{A total number openning hours during summer period}
#' \item{openning_hours}{A total number openning hours (winter + summer)}
#' \item{runway_area}{A total area of airport runways}
#' \item{terminal_area}{A total area of airport terminal(s)}
#' \item{parking_area}{A total area of airport parking area}
#' \item{island}{1 if an airpiort is located on an island; 0 otherwise}
#' \item{international}{1 if an airpiort is international; 0 otherwise}
#' \item{mixed_use}{1 if an airpiort is in mixed use; 0 otherwise}
#' \item{WLU}{A total volume of work load units (WLU) served by an airport}
#' \item{NearestCity}{A road network distance between an airport and its nearest city }
#' 
#' }
#' @source 
#' \describe{
#' \item{conference_paper}{"Airport efficiency and public investment in Greece" (2010) 
#' In Proceeding of the 2010 International Kuhmo-Nectar Conference on Transport Economics, University of Valencia, Spain.}
#' }
NULL


#' European airports statistical data (used for the RTE paper)
#' 
#' The \code{spfrontier} package includes the dataset \code{RTEpaper}, 
#' containing information about European airports infrastructure and traffic statistics in 2011 and used in 
#' Pavlyuk, D., 2016. Implication of spatial heterogeneity for airports’ efficiency estimation. 
#' Research in Transportation Economics 56, 15–24. https://doi.org/10.1016/j.retrec.2016.07.002.
#' 
#' 
#' @name RTEpaper
#' @rdname RTEpaper
#' @docType data
#' @format A dataframe with 39 observations on the following 24 variables.
#'
#' \describe{
#' 
#' \item{ID}{Airport identifier}
#' \item{ICAO_code}{Airport ICAO code}
#' \item{Country}{Airport country}
#' \item{AirportName}{Airport title}
#' \item{latitude}{Airport latitude}
#' \item{longitude}{Airport longitude}
#' \item{Year}{Year of observation}
#' \item{PAX}{A number of passengers served by an airport}
#' \item{cargo}{A total volume of cargo served by an airport}
#' \item{Population100km}{A number of inhabitants, living in 100 km around an airport}
#' \item{Island}{1 if an airport is located on an island; 0 otherwise}
#' \item{SouthIsland}{1 if an airport is located on a European South island; 0 otherwise}
#' \item{Routes}{Number of routes, served byt the airport}
#' \item{RunwayCount}{A number of airport runways}
#' \item{NUTS0}{NUTS0 region}
#' \item{NUTS1}{NUTS1 region}
#' \item{NUTS2}{NUTS2 region}
#' \item{NUTS3}{NUTS3 region}
#' \item{Tourists}{Number of tourists in airport's NUTS2 region}
#' \item{Area}{Area of airport's NUTS2 region}
#' \item{GDPpc}{Gross domestic product per capita in airport's NUTS3 region}
#' \item{hub}{Hub status of the airport, based on the number of served routes}
#' \item{international}{INternational status of the airport}
#' \item{ownership}{Ownership of the airport - public, minor private, major private, private}
#' \item{RunwayLength}{Total length of the airport's runways}
#' 
#' 
#' }
#' @source 
#' \describe{
#' \item{RTEpaper}{#' Pavlyuk, D., 2016. Implication of spatial heterogeneity for airports’ efficiency estimation. 
#' Research in Transportation Economics 56, 15–24. https://doi.org/10.1016/j.retrec.2016.07.002.}
#' }
NULL