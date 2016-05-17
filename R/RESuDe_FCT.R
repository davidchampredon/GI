################################################
###   Functions for Renewal Equation
###   with Susceptible Depletion (RESuDe)
###
###   Author: David Champredon
###
################################################


### Discrete generation interval distribution
###
gi.distrib <- function(type,mean,var,span,x){
	res <- NULL
	if(type == 'pois') res <- dpois(x = x,lambda = mean)
	if(type == 'nbinom') {
		if(var<mean) stop("For negBinom GI, we must have var>mean")
		res <- dnbinom(x = x,mu = mean, size = mean^2/(var-mean))
	}
	return(res)
}


### Simulation with RESuDe
### W a r n i n g: this implementation is deterministic
### because its use is to calculate generation
### interval distributions
###
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
