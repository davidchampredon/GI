###
###   Functions for Renewal Equation
###   with Susceptible Depletion (RESuDe)
###
###   Author: David Champredon
###
###


#' Internal function. Discrete generation interval distribution.
#' @param type String. Name of the distribution family for the GI. Implemented: \code{'pois'} or  \code{'nbinom'}
#' @param mean Numeric. Mean of the GI
#' @param var Numeric. Variance of the GI
#' @param span Numeric. Maximum value for the GI
#' @param x Numeric. Values where the GI must be calculated.
#' @return GI density distribution.
gi.distrib <- function(type,mean,var,span,x){
	res <- NULL
	if(type == 'pois') res <- dpois(x = x,lambda = mean)
	if(type == 'nbinom') {
		if(var<mean) stop("For negBinom GI, we must have var>mean")
		res <- dnbinom(x = x,mu = mean, size = mean^2/(var-mean))
	}
	return(res)
}


#' Internal function. Simulation with the Renewal Equation with Susceptible Depletion (RESuDe).
#' W a r n i n g: this implementation is deterministic because its use is to calculate generation interval distributions
#' @param pop_size Numeric. Population size.
#' @param I.init Numeric. Inial number of infectious individuals.
#' @param R0 Numeric. Basic reprocdution number.
#' @param alpha Numeric. Heterogeneous mixing parameter.
#' @param kappa Numeric. Intervention parameter.
#' @param GI_span Numeric. Maximum value for GI
#' @param GI_mean Numeric. Mean for GI
#' @param GI_var Numeric. Variance for GI
#' @param GI_type String. Type of distribution for the GI. \code{'pois'} or  \code{'nbinom'}\
#' @param horizon Numeric. Time horizon of the simulation.
RESuDe.simulate <- function(pop_size, 
							I.init,
							R0, 
							alpha, 
							kappa, 
							GI_span, 
							GI_mean, 
							GI_var,
							GI_type,
							horizon) {
	if(length(I.init)==1){
		# used to generate synthetic data
		I <- vector()
		S <- vector()
		I[1] <- I.init
		S[1] <- pop_size - I.init
	}
	if(length(I.init)>1){
	stop(paste("Parameter I.init must be a scalar, not a vector. Input value is:",I.init))
	}
	
	for(t in 2:(1+horizon)){
		z <- 0
		for(j in 1:min(GI_span,t-1)){
			GI_j <- gi.distrib(type = GI_type, 
							   mean = GI_mean, 
							   var = GI_var, 
							   span = GI_span,
							   x = j)
			z <- z + GI_j * I[t-j]
		}
		I.tmp <- (S[t-1]/ pop_size)^(1+alpha) * R0 * exp(-kappa*t) * z 
		I[t] <- min(I.tmp, S[t-1]) 
		S[t] <- max(0,S[t-1] - I[t])
	}
	return(list(S=S, I=I))
}
