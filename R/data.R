#' Backward generation intervals simulated by a SEmInR model
#'
#' A dataset containing the backward generation intervals simulated by a SEmInR model 
#' which is individual based and follows the Gillespie algorithm (with tau-leap) for the stochastic component.
#' The model is wrapped in a R package and can be downloaded here: \code{https://github.com/davidchampredon/seminribm}.
#' 
#' The parameter values used for this data set are:
#' horizon <- 300
#' popSize <- 5e3
#' 
#' initInfectious  <- 2
#' R0              <- 3.0
#' latent_mean     <- 2
#' infectious_mean <- 4
#' nE              <- 6
#' nI              <- 6
#' calc_WIW_Re     <- FALSE
#' doExact         <- FALSE
#' timeStepTauLeap <- 0.1
#' rnd_seed        <- 1234
#'
#' @format A data frame with 4749 rows (=infection events) and 3 variables:
#' \describe{
#'   \item{at}{Disease acquisition time.}
#'   \item{rt}{Rounded value of \code{at}.}
#'   \item{b}{Value of the backward generation interval for this infection event.}
#' }
#' @source \url{https://github.com/davidchampredon/seminribm}
"sim_gi"