### 
###   GENERATION INTERVAL FUNCTIONS
###
###   Author: David Champredon
###


#' Internal function. Calculate intrinsic generation interval distribution.
#' @param tau Integer. The _index_ of the simulation time (not the time value.)
#' @param df Dataframe of the epidemic simulation time series.
#' @return Density of the intrinsic genration interval distribution.
#' 
GI.intrinsic <- function(tau, df){
	res <- 0
	N  <- length(df$time)
	dt <- df$time[2]-df$time[1]
	JJ <- df$Jall[1:N]
	integ  <- sum(JJ)*dt
	res <- df$Jall[tau]/integ
	return(res)
}

#' Internal function. Calculate the _forward_ generation interval distribution.
#' @param tau Integer. The _index_ of the simulation time (not the time value.)
#' @param s Integer. The _index_ of the simulation time (not the time value.)
#' @param df Dataframe of the epidemic simulation time series.
#' @return Density of the forward genration interval distribution.
#' 
GI.fwd <- function(tau, s, df){
	
	# Warning: 'tau' and 's' is the _index_ of the simulation time 
	# (not the value of the time itself!)
	res <- 0
	N  <- length(df$time)
	dt <- df$time[2]-df$time[1]
	
	JJ <- df$Jall[1:(N-s)]
	SS <- df$S[(s+1):N]
	integ  <- sum(JJ*SS)*dt
	res <- df$Jall[tau]*df$S[tau+s]/integ
	return(res)
}

#' Internal function. Calculate the _backward_ generation interval distribution.
#' @param tau Integer. The _index_ of the simulation time (not the time value.)
#' @param t Integer. The _index_ of the simulation time (not the time value.)
#' @param df Dataframe of the epidemic simulation time series.
#' @return Density of the backward genration interval distribution.
#' 
GI.bck <- function(tau, t, df){
	
	# Warning: 'tau' and 't' is the _index_ of the simulation time 
	# (not the value of the time itself!)
	
	res <- 0
	
	if(tau<t){
		N <- length(df$time)
		dt <- df$time[2]-df$time[1]
		
		JJ <- df$Jall[1:(t-1)]
		II <- df$inc[(t-1):1]
		integ <- sum(JJ*II)*dt
		res <- df$Jall[tau]*df$inc[t-tau]/integ
	}
	return(res)
}


#' Internal function. Calculate the _forward_ generation interval distribution with discrete times.
#' @param time.s Integer. Calendar time.
#' @param g Numerical vector. Intrinsic generation interval distribution.
#' @param GI_span Maximum value for the intrinsic GI.
#' @param S Numerical vector. Susceptible time series.
#' @return Density of the forward genration interval distribution.
#' 
GI.fwd.discrete <- function(time.s, g, GI_span, S) {
	gfwd <- vector()
	for (tau in 1:GI_span) {
		tmp <- 0
		kmax <- min(length(g), length(S)-time.s+1)
		for (k in 1:length(g))
			tmp <- tmp + g[k] * S[time.s + k]
		gfwd[tau] <- g[tau] * S[time.s + tau] / tmp
	}
	return(gfwd)
}

#' Internal function. Calculate the _backward_ generation interval distribution with discrete times.
#' @param time.t Integer. Calendar time.
#' @param g Numerical vector. Intrinsic generation interval distribution.
#' @param GI_span Maximum value for the intrinsic GI.
#' @param I Numerical vector. Prevalence time series.
#' @return Density of the backward genration interval distribution.
#' 
GI.bck.discrete <- function(time.t, g, GI_span, I) {
	gbck <- vector()
	for (tau in 1:min(GI_span, time.t - 1)) {
		tmp <- 0
		kmax <- min(length(g), time.t - 1)
		for (k in 1:kmax)
			tmp <- tmp + g[k] * I[time.t - k]
		gbck[tau] <- g[tau] * I[time.t - tau] / tmp
	}
	return(gbck)
}

#' @title Generation interval distribution for the RESuDe epidemic model.
#' @description  Calculates the forward and backward generation interval distributions for the Renwal Equation with Susceptible Depletion (RESuDe) epidemic model, defined as:
#' i(t) = R0*(S(t)/N)^{alpha} * exp(-kappa*t) sum_{k=1}^{GI_span} g(k) * i(t-k)
#' S(t) = max(0; S(t-1) - i(t) )
#' with S the number of susceptible individuals, i the incidence, N the constant population size. The RESuDe model is specified directly with the intrinsic generation interval distribution, g, but it is not straightforward to calculate the backward and forward generation interval distributions. This function numericaly solves the backward and forward generation interval distributions (see [1] for theoretical framework).
#' @param cal.times.fwdbck Numeric vector. Times where the backward and forward generation interval distributions will be calculated
#' @param R0 Numeric. The basic reproduction number
#' @param alpha Numeric. Mixing heterogeneity parameter
#' @param kappa Numeric. Exponential rate decreasing the force of infection (for example because of interventions)
#' @param GI_span Numeric.  The support for the intrinsic generation interval distribution (the maximum value of the generation interval)
#' @param GI_mean Numeric. Mean of the _intrinsic_ generation interval
#' @param GI_var Numeric. Variance of the _intrinsic_ generation interval
#' @param GI_type String. Distribution type for the generation interval. \code{pois} for Poisson, \code{nbinom} for negative binomial.
#' @param horizon Numeric.  Time until which the simultion is run. User must make sure it is beyond the end of the epidemic, else the calculation of generation interval distributions may be wrong!
#' @param pop_size Integer. Constant population size. Default value at 1E4
#' @param I.init Integer. Initial number of infectious individuals. Default value at 1
#' @return Returns a named list.
#' \itemize{
#' \item{intrinsic:} {Intrinsic generation interval distribution. List composed of two elements: tsi, the time since infection vector and density the associated vector of densities. }
#' \item{fwd:} {Forward generation interval distribution. List composed of cal.times.fwdbck elements, each one representing the densities of the forward generation interval distribution at the requested calendar times.}
#' \item{bck:} Backward generation interval distribution. List composed of cal.times.fwdbck elements, each one representing the densities of the backward generation interval distribution at the requested calendar times. 
#' \item{fwd.mean:} Mean of the forward generation interval distribution. List composed of cal.times.fwdbck elements, each one representing the mean of the forward generation interval distribution at the requested calendar times.
#' \item{bck.mean:} Mean of the backward generation interval distribution. List composed of cal.times.fwdbck elements, each one representing the mean of the backward generation interval distribution at the requested calendar times.
#' \item{incidence:} Incidence modeled from the RESuDe model.
#' \item{susceptible:} Time series of the number of susceptible individuals modeled from the RESuDe model.
#' }
#' @export
#' @references [1] Champredon D, Dushoff J. Intrinsic and realized generation intervals in infectious-disease transmission. Proceedings of the Royal Society B: Biological Sciences 2015; 282: 20152026.
#'
GI.resude <- function(cal.times.fwdbck,
					  R0,
					  alpha,
					  kappa,
					  GI_span,
					  GI_mean,
					  GI_var,
					  GI_type,
					  horizon,
					  pop_size = 1E4,
					  I.init = 1){
	
	sim <- RESuDe.simulate(pop_size,
						   I.init,
						   R0,
						   alpha,
						   kappa,
						   GI_span,
						   GI_mean,
						   GI_var,
						   GI_type,
						   horizon)
	
	xx <- 1:GI_span
	g <- gi.distrib(
			type = GI_type,
			mean = GI_mean,
			var = GI_var,
			span = GI_span,
			x = xx
		)
	S <- sim$S
	I <- sim$I
	
	gfwd <- list()
	gbck <- list()
	fbar <- vector()
	bbar <- vector()
	cnt <- 1
	for(cal.t in cal.times.fwdbck){
		gfwd[[cnt]] <- GI.fwd.discrete(time.s=cal.t, g, GI_span, S)
		gbck[[cnt]] <- GI.bck.discrete(time.t=cal.t, g, GI_span, I)
		nf <- 1:length(gfwd[[cnt]])
		nb <- 1:length(gbck[[cnt]])
		fbar[cnt] <- sum(nf*gfwd[[cnt]])
		bbar[cnt] <- sum(nb*gbck[[cnt]])
		cnt <- cnt + 1
	}
	
	return(list(intrinsic = list(tsi = (1:length(g)), density =g),
				fwd       = list(tsi = (1:length(g)), density = gfwd),
				bck       = list(tsi = (1:length(g)), density = gbck),
				fwd.mean  = fbar,
				bck.mean  = bbar,
				incidence = I,
				susceptible = S)
	)
}

#' Calculates the intrinsic, forward and backward generation interval distributions for SEmInR epidemic model. 
#' The SEmInR model is an epidemiological model that compartments individuals by stages of infections: S for susceptible, E for exposed (infected but not yet infectious), I infectious and R for individuals removed from the model (following recovery and immunity, deaths, etc.). In order to have realistic residency distributions in compartments E and I, these compartments are artificially sub-divided into m for E, and n for I. 
#' It is not straightforward to calculate the generation interval distributions for a SEmInR model. This function numericaly solves the three generation interval distributions (intrinsic, backward and forward ones, see [1] for theoretical framework).
#' @param cal.times.fwdbck Numeric vector. Times where the backward and forward generation interval distributions will be calculated
#' @param R0 Numeric. The basic reproduction number
#' @param latent_mean Numeric. The mean duration of the latent period (in 'E' compartments)
#' @param infectious_mean Numeric. The mean duration of the infectious period (in 'I' compartments)
#' @param nE Integer Number of E compatments
#' @param nI Integer. Number of I compatments
#' @param horizon Numeric.  Time until which the simultion is run. User must make sure it is beyond the end of the epidemic, else the calculation of generation interval distributions may be wrong!
#' @param pop_size Integer. Constant population size. Default value at 1E4
#' @param I.init Integer. Initial proportion of infectious individuals. Default value at 1
#' @return Returns a named list.
#' \itemize{
#' \item{intrinsic:} {Intrinsic generation interval distribution. List composed of two elements: tsi, the time since infection vector and density the associated vector of densities. }
#' \item{fwd:} {Forward generation interval distribution. List composed of cal.times.fwdbck elements, each one representing the densities of the forward generation interval distribution at the requested calendar times.}
#' \item{bck:} Backward generation interval distribution. List composed of cal.times.fwdbck elements, each one representing the densities of the backward generation interval distribution at the requested calendar times. 
#' \item{fwd.mean:} Mean of the forward generation interval distribution. List composed of cal.times.fwdbck elements, each one representing the mean of the forward generation interval distribution at the requested calendar times.
#' \item{bck.mean:} Mean of the backward generation interval distribution. List composed of cal.times.fwdbck elements, each one representing the mean of the backward generation interval distribution at the requested calendar times.
#' \item{incidence:} Incidence modeled from the SEmInR model.
#' \item{susceptible:} Time series of the number of susceptible individuals modeled from the SEmInR model.
#' }
#' @references [1] Champredon D, Dushoff J. Intrinsic and realized generation intervals in infectious-disease transmission. Proceedings of the Royal Society B: Biological Sciences 2015; 282: 20152026.
#' @importFrom deSolve lsoda
#' @export
#' 
GI.seminr <- function(latent_mean,
					  infectious_mean,
					  R0,
					  nE, 
					  nI,
					  cal.times.fwdbck,
					  horizon,
					  dt = 0.1,
					  I.init = 1E-5
){
	# Generation interval distributions have to
	# be solved numerically for SEmInR models
	# (no analytical formula). Hence, first simulate
	# the epidemic, then calculate generation
	# interval distributions.
	
	############################
	#### SEmInR Simulations ####
	############################
	
	# SEmInR Parameters
	tvec  <- seq(0, horizon, by = dt)
	sigma <- 1/latent_mean
	gamma <- 1/infectious_mean
	beta  <- R0 * gamma
	f     <- 0.0
	mu    <- 0.0
	
	# Initial conditions
	S.init <- 1 - I.init
	E.init <- 0
	
	params.SEmInR <- c(mu=mu, 
					   beta=beta,
					   sigma=sigma,
					   gamma=gamma,
					   f=f,
					   nE=nE, nI=nI)
	
	# Inital conditions
	# W A R N I N G : initial infectious in I[1] compartment (not E[1])
	EIvec <- c(E= rep(0,ifelse(nE==Inf,1,nE)),
			   I= c(I.init,rep(0,ifelse(nI==Inf,0,nI-1))) )
	
	YJvec <- c(Y=c(1,rep(0,ifelse(nE==Inf,1,nE-1))),
			   J=rep(0,ifelse(nI==Inf,0,nI) ))
	
	inits.SEmInR <- c(S=1-I.init, EIvec, R=0, YJvec, Z=I.init)
	
	# Solutions:
	SEmInR <- as.data.frame(lsoda(y     = inits.SEmInR, 
								  times = tvec, 
								  func  = SEmInR.FX, 
								  parms = params.SEmInR))
	SEmInR <- calc.Iall(dat = SEmInR)
	SEmInR <- calc.Jall(dat = SEmInR)
	
	# incidence (from solved cumulative incidence)
	SEmInR$inc <- c(SEmInR$Z[1],diff(SEmInR$Z))
	
	##########################
	#### Generation times ####
	##########################
	
	NT <- nrow(SEmInR)
	tsim <- SEmInR$time
	
	# mean forward & backward generation interval
	tt <- vector()
	f.bar <- vector()
	g.bar <- vector()
	
	# Calculate the indices where simulation times are 
	# close to requested calendar times for fwd & bck GI distributions:
	cal.t.idx <- vector()
	for(k in 1:length(cal.times.fwdbck)) {
		cal.t.idx[k] <- which(abs(tsim-cal.times.fwdbck[k])<dt/1.99)[1]
	}
	
	# calculate gi.intrinsic
	gi.intrinsic <- vector()
	for(tau in 1:NT) { gi.intrinsic[tau] <- GI.intrinsic(tau,SEmInR) }
	
	# theoretical mean: 
	theo.mean.gii <- latent_mean + infectious_mean*(nI+1)/2/nI
	# numerical mean:
	mean.gii <- sum(SEmInR$time*gi.intrinsic*dt)
	# variance:
	var.gii <- sum(SEmInR$time^2*gi.intrinsic*dt) - mean.gii^2
	
	### Forward and Backward generation interval distributions
	cnt <- 1
	gi.fwd.list <- list()
	gi.bck.list <- list()
	gi.fwd.t.list <- list()
	gi.bck.t.list <- list()
	
	for(s in cal.t.idx){
		
		# calculate gi.fwd:
		gi.fwd <- vector()
		for(tau in 1:(NT-s)) { gi.fwd[tau] <- GI.fwd(tau, s, SEmInR) }
		# calculate expectation of gi.fwd:
		tt[cnt] <- SEmInR$time[s]
		tvec <- SEmInR$time[c(1:(NT-s))]
		f.bar[cnt] <- sum(gi.fwd*tvec*dt)
		
		# calculate gi.bck:
		gi.bck <- vector()
		for(tau in 1:(s-1)){ gi.bck[tau] <- GI.bck(tau, s, SEmInR) }  
		# calculate expectation of gi.bck:
		tvec.bck <- SEmInR$time[c(1:(s-1))]
		g.bar[cnt] <- sum(gi.bck*tvec.bck*dt)
		
		# Store distribution at each calendar time requested:
		gi.fwd.list[[cnt]]   <- gi.fwd
		gi.fwd.t.list[[cnt]] <- tvec[1:length(gi.fwd)]
		gi.bck.list[[cnt]]   <- gi.bck[-length(gi.bck)] # remove last value bc numerical instability
		gi.bck.t.list[[cnt]] <- tvec.bck[1:(length(gi.bck)-1)]
		
		cnt <- cnt+1
	}
	return(list(intrinsic = list(tsi = tsim,          density = gi.intrinsic), # tsi=time since infection
				fwd       = list(tsi = gi.fwd.t.list, density = gi.fwd.list),
				bck       = list(tsi = gi.bck.t.list, density = gi.bck.list),
				fwd.mean  = f.bar,
				bck.mean  = g.bar,
				incidence = SEmInR$inc,
				cal.times = tsim))
}

