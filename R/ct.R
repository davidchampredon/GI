### 
###   Fit to tontact tracing data
###

#' Negative Log-Likelihood given R0 and mean generation interval
#' assuming Poisson observation error from contact tracing data
#' 
#' @param t.obs Numeric vector. Time when the backward generation intervals were observed
#' @param gi.obs Numeric vector. Backward generation intervals observed.
#' @param model.epi String. Epidemic model used. Choice between \code{seminr} or \code{resude}.
#' @param fxd.prm List. Parameters of \code{model.epi} that are fixed.
#' 
nllk <- function(R0, 
                 gimean,
                 t.obs, 
                 gi.obs, 
                 model.epi, 
                 fxd.prm) {
    
    if(model.epi=='seminr'){
        
        im <- 2*(gimean - fxd.prm[['latent_mean']])  # gi.mean ~ latent + infectious/2
        if(im <= 0){
            warning('GI and latent period lengths not consistent!')
            im <- 1
        }
        G <- GI.seminr(latent_mean = fxd.prm[['latent_mean']],
                       infectious_mean = im, 
                       R0 = R0, 
                       nE = fxd.prm[['nE']], 
                       nI = fxd.prm[['nI']],
                       cal.times.fwdbck = t.obs,
                       horizon = fxd.prm[['horizon']], 
                       dt = fxd.prm[['dt']])
    }
    if(model.epi=='resude'){
        G <- GI.resude(cal.times.fwdbck = t.obs,
                       R0 = R0,
                       alpha = fxd.prm[['alpha']], 
                       kappa = fxd.prm[['kappa']], 
                       GI_span = fxd.prm[['GI_span']], 
                       GI_mean = gimean, 
                       GI_var = fxd.prm[['GI_var']], 
                       GI_type = fxd.prm[['GI_type']],
                       horizon = fxd.prm[['horizon']])
    }
    # Extract the mean bckwd GI :
    b <- G$bck.mean    
    b[b==0] <- 1e-6
    
    # Calculate log likelihood:
    z <- sum(dpois(x = gi.obs, lambda = b, log = TRUE))
    return(-z)
}


#' Fit R0 and mean generation interval from contact tracing data
#' 
#' @param t.obs Numeric vector. Time when the backward generation intervals were observed
#' @param gi.obs Numeric vector. Backward generation intervals observed.
#' @param model.epi String. Epidemic model used. Choice between \code{seminr} or \code{resude}.
#' @param fxd.prm List. Parameters of \code{model.epi} that are fixed.
#' @param R0.rng Numeric vector. Values of R0 explored to calculate the likelihood surface.
#' @param gimean.rng Numeric vector. Values of mean intrinsic GI explored to calculate the likelihood surface.
#' @param CI Numeric. Confidence interval level. Default = 0.95.
#' @param do.plot Boolean. Fitting diagnostic plots. Default = TRUE.
#' @export
gi_ct_fit <- function(t.obs, 
                      gi.obs, 
                      model.epi, 
                      fxd.prm,
                      R0.rng, 
                      gimean.rng,
                      CI = 0.95,
                      do.plot = FALSE ) {
    
    t1 <- as.numeric(Sys.time())
    
    M <- outer(X = R0.rng, 
               Y = gimean.rng, 
               FUN = Vectorize(nllk, list("R0","gimean")), 
               t.obs = t.obs, 
               gi.obs = gi.obs, 
               model.epi = model.epi, 
               fxd.prm = fxd.prm)
    
    idx <- which(M == min(M, na.rm = TRUE), 
                 arr.ind = TRUE)
    
    # Retrieve neg log likelihood min
    ll.min <- M[idx]
    
    # Calculate the level on the likelihoof function
    conf.cutoff <- ll.min + qchisq(CI,2)/2
    
    # Estimate the confidence interval from the contour lines
    cc <- contourLines(R0.rng,
                       gimean.rng,
                       M,
                       level = conf.cutoff)
    x.tmp <- list()
    y.tmp <- list()
    
    for(i in seq_along(cc)){
        x.tmp[[i]] <- range(cc[[i]]$x)
        y.tmp[[i]] <- range(cc[[i]]$y)
    }
    R0.ci     <- range(unlist(x.tmp))
    gimean.ci <- range(unlist(y.tmp))
    
    # Likelihood surface:
    if(do.plot){
        title <- paste0('Neg. Log-Likelihood Surface (',model.epi,')\n with ',CI*100,'%CI')
        ylab <- 'Mean intrinsic GI'
        if(model.epi=='seminr') ylab <- 'Mean infectious period'
        
        contour(x = R0.rng, 
                y = gimean.rng, 
                z = M, 
                nlevels = 20,
                col='grey',
                main = title,
                xlab = 'R0', 
                ylab = ylab,
                las = 1)
        
        # Best point estimate:
        col.best <- 'red'
        R0.best     <- R0.rng[idx[1]]
        gimean.best <- gimean.rng[idx[2]]
        
        points(x = R0.best, 
               y = gimean.best,
               pch=16, cex=3, col=col.best)
        segments(x0 = R0.best, 
                 y0 = 0,
                 x1 = R0.best,
                 y1 = gimean.best,
                 col = col.best, lty=2)
        segments(x0 = 0,
                 y0 = gimean.best,
                 x1 = R0.best,
                 y1 = gimean.best,
                 col = col.best, lty=2)
        
        # Draw the confidence contour
        contour(R0.rng,
                gimean.rng,
                M,
                level = conf.cutoff,
                labels="",
                col="red",
                lwd=3, 
                lty=1,
                add = TRUE)
        
        # Plot the fitted GI means:
        
        if(model.epi=='seminr'){
            im <- 2*(gimean.best - fxd.prm[['latent_mean']])  # gi.mean ~ latent + infectious/2
            im.lo <- 2*(gimean.ci[1] - fxd.prm[['latent_mean']])  # gi.mean ~ latent + infectious/2
            im.hi <- 2*(gimean.ci[2] - fxd.prm[['latent_mean']])  # gi.mean ~ latent + infectious/2
            if(im <= 0) im <- 1
            if(im.lo <= 0) im.lo <- 1
            if(im.hi <= 0) im.hi <- 1
            Gfit <- GI.seminr(latent_mean = fxd.prm[['latent_mean']],
                              infectious_mean = im, 
                              R0 = R0.best, 
                              nE = fxd.prm[['nE']], 
                              nI = fxd.prm[['nI']],
                              cal.times.fwdbck = t.obs,
                              horizon = fxd.prm[['horizon']], 
                              dt = fxd.prm[['dt']])
            
            Gfit.lo <- GI.seminr(latent_mean = fxd.prm[['latent_mean']],
                              infectious_mean = im.lo, 
                              R0 = R0.ci[1], 
                              nE = fxd.prm[['nE']], 
                              nI = fxd.prm[['nI']],
                              cal.times.fwdbck = t.obs,
                              horizon = fxd.prm[['horizon']], 
                              dt = fxd.prm[['dt']])
            Gfit.hi <- GI.seminr(latent_mean = fxd.prm[['latent_mean']],
                                 infectious_mean = im.hi, 
                                 R0 = R0.ci[2], 
                                 nE = fxd.prm[['nE']], 
                                 nI = fxd.prm[['nI']],
                                 cal.times.fwdbck = t.obs,
                                 horizon = fxd.prm[['horizon']], 
                                 dt = fxd.prm[['dt']])
            
        }
        if(model.epi=='resude'){
            Gfit <- GI.resude(cal.times.fwdbck = t.obs,
                              R0 = R0.best,
                              alpha = fxd.prm[['alpha']], 
                              kappa = fxd.prm[['kappa']], 
                              GI_span = fxd.prm[['GI_span']], 
                              GI_mean = gimean.best, 
                              GI_var = fxd.prm[['GI_var']], 
                              GI_type = fxd.prm[['GI_type']],
                              horizon = fxd.prm[['horizon']])
            Gfit.lo <- GI.resude(cal.times.fwdbck = t.obs,
                              R0 = R0.ci[1],
                              alpha = fxd.prm[['alpha']], 
                              kappa = fxd.prm[['kappa']], 
                              GI_span = fxd.prm[['GI_span']], 
                              GI_mean = gimean.ci[1], 
                              GI_var = fxd.prm[['GI_var']], 
                              GI_type = fxd.prm[['GI_type']],
                              horizon = fxd.prm[['horizon']])
            Gfit.hi <- GI.resude(cal.times.fwdbck = t.obs,
                                 R0 = R0.ci[2],
                                 alpha = fxd.prm[['alpha']], 
                                 kappa = fxd.prm[['kappa']], 
                                 GI_span = fxd.prm[['GI_span']], 
                                 GI_mean = gimean.ci[2], 
                                 GI_var = fxd.prm[['GI_var']], 
                                 GI_type = fxd.prm[['GI_type']],
                                 horizon = fxd.prm[['horizon']])
        }
        
        gbck.fit <- Gfit$bck.mean
        gbck.fit.lo <- Gfit.lo$bck.mean
        gbck.fit.hi <- Gfit.hi$bck.mean
        
        plot(x = t.obs, 
             y = gbck.fit, 
             ylim = range(gbck.fit,gi.obs),
             typ='o', pch=16, lwd=3,
             col = col.best,
             las = 1, 
             xlab = 'calendar time',
             ylab='Mean Backward GI',
             main = paste(model.epi, 'model backward GI\nfitted to contact tracing data'))
        lines(x=t.obs, y=gbck.fit.lo, lty=2, col=col.best, lwd=2)
        lines(x=t.obs, y=gbck.fit.hi, lty=2, col=col.best, lwd=2)
        points(t.obs, gi.obs, pch=1, col='black',lwd=1.5,cex=1)
        grid()
        legend('topleft', legend = c('data','model fit',paste(CI*100,'%CI')),
               col=c('black',col.best,col.best),pch=c(1,16,NA),lwd=c(NA,3,2),
               pt.cex = c(1,1,NA), pt.lwd = c(1.5,1,NA), lty=c(1,1,2))
    }
    t2 <- as.numeric(Sys.time())
    dt <- round( (t2-t1)/60, 1)
    msg <- paste('Fit to contact tracing data done in',dt,'minute(s).')
    message(msg)
    return(list(R0.best = R0.best,
                gimean.best = gimean.best,
                R0.ci = R0.ci,
                gimean.ci = gimean.ci))
    
}