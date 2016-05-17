
GI.intrinsic <- function(tau, df){
	
	# Warning: 'tau' is the _index_ of the simulation time 
	# (not the value of the time itself!)
	res <- 0
	N  <- length(df$time)
	dt <- df$time[2]-df$time[1]
	JJ <- df$Jall[1:N]
	integ  <- sum(JJ)*dt
	res <- df$Jall[tau]/integ
	return(res)
}

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
	# the epidemic, then caluclate generation
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
	SEmInR <- calc.Yall(dat = SEmInR)
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

