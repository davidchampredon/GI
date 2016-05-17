\name{GI.seminr}
\alias{GI.seminr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation interval distributions for SEmInR epidemic model
}
\description{
Calculates the intrinsic, forward and backward generation interval distributions for SEmInR epidemic model. 
The SEmInR model is an epidemiological model that compartments individuals by stages of infections: S for susceptible, E for exposed (infected but not yet infectious), I infectious and R for individuals removed from the model (following recovery and immunity, deaths, etc.). In order to have realistic residency distributions in compartments E and I, these compartments are artificially sub-divided into m for E, and n for I. 
It is not straightforward to calculate the generation interval distributions for a SEmInR model. This function numericaly solves the three generation interval distributions (intrinsic, backward and forward ones, see [1] for theoretical framework).
}
\usage{
GI.seminr(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{latent_mean}{
  : The mean duration of the latent period (in 'E' compartments)
}
  \item{infectious_mean}{
  : The mean duration of the infectious period (in 'I' compartments)
}
  \item{R0}{
  : The basic reproduction number
}
  \item{nE}{
  : Number of E compatments
}
  \item{nI}{
  : Number of I compatments
}
  \item{cal.times.fwdbck}{
  : Calendar times where the backward and forward generation interval distributions will be calculated
}
  \item{dt}{ 
  : Time step for ODE solver. Defaults to 0.1
}
\item{I.init}{ 
  : Initial proportion of infectious individual to the total population. Defaults to 1E-5
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
: Incidence modeled from the SEmInR model.
}
\item{cal.times}{
: Vector of (calendar) times used in the numerical schemes.
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
}
\examples{
library(GI)


# SEmInR Parameters
latent_mean     <- 5
infectious_mean <- 5
R0 <- 4
nE <- 3
nI <- 3
horizon <- 150
dt      <- 0.1
I.init  <- 1E-5

cal.times.fwdbck <- seq(2, horizon-4, by=2)

g <- GI.seminr(latent_mean, infectious_mean, R0, nE,nI, cal.times.fwdbck,horizon,dt,I.init)
g.int <- g[['intrinsic']]
g.fwd <- g[['fwd']]
g.bck <- g[['bck']]
t     <- g[['cal.times']]
inc   <- g[['incidence']]

# PLOTS 
par(mfrow=c(2,3))

plot(t,inc, typ='l', 
	 main = "SEmInR simulated incidence",
	 xlab = 'calendar time')

plot(g.int$tsi, g.int$density, 
	 typ='l', col='red',
	 xlab = 'time since infection', ylab ='',
	 main='Intrinsic generation interval distribution')

plot(x=g.fwd$tsi[[1]], y = g.fwd$density[[1]],
	 main = 'Forward generation interval distribution',
	 xlab = 'time since infection', ylab ='',
	 typ='l',ylim=c(0,0.4))
for(i in 1:length(g.fwd$tsi)) lines(x=g.fwd$tsi[[i]],g.fwd$density[[i]])

plot(x=g.bck$tsi[[length(g.bck$tsi)]], g.bck$density[[length(g.bck$tsi)]],
	 main = 'Backward generation interval distribution',
	 xlab = 'time since infection', ylab ='',
	 typ='l',ylim=c(0,0.4))
for(i in 1:length(g.bck$tsi)) lines(x=g.bck$tsi[[i]],g.bck$density[[i]])

fbar <- g[['fwd.mean']]
plot(cal.times.fwdbck,fbar,typ='o', 
	 main = 'Mean of the forward generation interval distribution',
	 xlab = 'calendar time', ylab ='')

bbar <- g[['bck.mean']]
plot(cal.times.fwdbck,bbar,typ='o',
	 main = 'Mean of the backward generation interval distribution',
	 xlab = 'calendar time', ylab ='')

}