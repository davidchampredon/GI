\name{GI.resude}
\alias{GI.resude}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation interval distributions for RESuDe epidemic model
}
\description{
Calculates the forward and backward generation interval distributions for the Renwal Equation with Susceptible Depletion (RESuDe) epidemic model, defined as:

\deqn{ i(t) = R0*(S(t)/N)^{alpha} * exp(-kappa*t) sum_{k=1}^{GI_span} g(k) * i(t-k)  }
\deqn{ S(t) = max(0; S(t-1) - i(t) ) }
with S the number of susceptible individuals, i the incidence, N the constant population size. The RESuDe model is specified directly with the intrinsic generation interval distribution, g, but it is not straightforward to calculate the backward and forward generation interval distributions. This function numericaly solves the backward and forward generation interval distributions (see [1] for theoretical framework).
}
\usage{
GI.resude(cal.times.fwdbck, R0, alpha, kappa, GI_span, GI_mean, GI_var, GI_type, horizon, pop_size = 1E4, I.init = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{cal.times.fwdbck}{
  : Calendar times where the backward and forward generation interval distributions will be calculated
}
\item{R0}{
  : The basic reproduction number
}
\item{alpha}{
  : Mixing heterogeneity parameter
}
\item{kappa}{
  : Exponential rate decreasing the force of infection (for example becaus of interventions)
}
\item{GI_span}{
  : The support for the intrinsic generation interval distribution (the maximum value of the generation interval)
}
\item{GI_mean}{
  : Mean of the intrinsic generation interval
}
\item{GI_var}{
  : Variance of the intrinsic generation interval
}
\item{GI_type}{
  : Distribution type for the generation interval. \code{pois} for Poisson, \code{nbinom} for negative binomial.
}
\item{horizon}{
  : Time until which the simultion is run. \emph{User must make sure it is beyond the end of the epidemic, else the calculation of generation interval distributions may be wrong!}
}
\item{pop_size}{
  : Constant population size. Default value at 1E4
}
\item{I.init}{
  : Initial number of infectious individuals. Default value at 1
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

Returns a list:
\item{intrinsic}{
: Intrinsic generation interval distribution. List composed of two elements: \code{tsi}, the time since infection vector and \code{density} the associated vector of densities.
}
\item{fwd}{
: Forward generation interval distribution. List composed of \code{cal.times.fwdbck} elements, each one representing the densities of the forward generation interval distribution at the requested calendar times.
}
\item{bck}{
: Backward generation interval distribution. List composed of \code{cal.times.fwdbck} elements, each one representing the densities of the backward generation interval distribution at the requested calendar times.
}
\item{fwd.mean}{
: Mean of the forward generation interval distribution. List composed of \code{cal.times.fwdbck} elements, each one representing the mean of the forward generation interval distribution at the requested calendar times.
}
\item{bck.mean}{
: Mean of the backward generation interval distribution. List composed of \code{cal.times.fwdbck} elements, each one representing the mean of the backward generation interval distribution at the requested calendar times.
}
\item{incidence}{
: Incidence modeled from the RESuDe model.
}
\item{susceptible}{
: Time series of the number of susceptible individuals modeled from the RESuDe model.
}
}

\references{
%% ~put references to the literature/web site here ~
[1] Champredon D, Dushoff J. Intrinsic and realized generation intervals in infectious-disease transmission. Proceedings of the Royal Society B: Biological Sciences 2015; 282: 20152026.
}
\author{
David Champredon
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{GI.seminr}}
}
\examples{
library(GI)

### Parameters for RESuDe model
R0 <- 4
horizon <- 100
I.init <- 1
pop_size <- 1e4
alpha <- 0
kappa <- 0
GI_type <- 'pois'  # 'pois' or 'nbinom'
GI_span <- 20
GI_mean <- 8
GI_var <- GI_mean + 2  # used only if GI_type='nbinom'

# Calendar times where the forward and 
# backward generation interval distributions
# will be calculated:
cal.times.fwdbck <- seq(GI_span-2, horizon - 4, by = 1)

# Calculate generation interval distributions:
G <- GI.resude(cal.times.fwdbck,
				R0,
				alpha,
				kappa,
				GI_span,
				GI_mean,
				GI_var,
				GI_type,
				horizon,
				pop_size = pop_size,
				I.init = I.init)

# Retrieve results:
gi <- G$intrinsic
gbck <- G$bck
gfwd <- G$fwd


###  PLOTS  ###

par(mfrow=c(3,2)) 

# Epidemic time series
plot(G$susceptible,
	 xlab = 'calendar time', ylab='', 
	 main='Number of susceptibles',
	 typ='l')
grid()
plot(G$incidence,typ='l',
	 xlab = 'calendar time', ylab='', 
	 main = 'Incidence')
grid()

# GI Distributions:
plot(x=gi$tsi, y=gi$density,
	 main = 'Backward GI distributions\n at requested calendar times',
	 xlab = 'time since infection', ylab='',yaxt='n',
	 ylim = c(0,1.2*max(gi$density)),
	 pch=16)
for(i in 1:length(cal.times.fwdbck)){
	gb <- gbck$density[[i]]
	lines(x = gbck$tsi[1:length(gb)], y = gb, typ='l', col = rgb(0,0,1,0.2))
}
legend(x='topright',lwd = c(1,0),legend = c('backward','intrinsic'),pch=c(NA,16))

plot(x=gi$tsi, y=gi$density,
	 main = 'Forward GI distributions\n at requested calendar times',
	 xlab = 'time since infection', ylab='',yaxt='n',
	 ylim = c(0,1.2*max(gi$density)),
	 pch=16)
for(i in 1:length(cal.times.fwdbck)){
	lines(x = gfwd$tsi, y = gfwd$density[[i]], typ='l', col=rgb(0,.5,0,0.1))
}
legend('topright',lwd = c(1,0),legend = c('forward','intrinsic'),pch=c(NA,16))


plot(x=cal.times.fwdbck, y=G$fwd.mean, 
	 xlim = c(0,horizon),
	 main = 'Mean of the forward GI distribution\n at requested calendar times',
	 xlab = 'calendar time', ylab='', 
	 las = 1,
	 typ='o')
abline(h=GI_mean, lty=2)
grid()
plot(x=cal.times.fwdbck, y=G$bck.mean, 
	 xlim = c(0,horizon),
	 main = 'Mean of the forward GI distribution\n at requested calendar times',
	 xlab = 'calendar time', ylab='', 
	 las = 1,
	 typ='o')
abline(h=GI_mean, lty=2)
grid()

}