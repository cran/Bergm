\name{bergm}
\Rdversion{1.2}
\alias{bergm}
\title{
Bayesian exponential random graph models
}
\description{
Function to fit Bayesian exponential random graphs models
using the exchange algorithm with either single-site update
or population MCMC with snooker update.
}
\usage{
bergm(model, burn.in = 0, main.iter = 5000, aux.iter = 1000, 
sdprop = NULL, sdprior = NULL, mprior = NULL, 
popMCMC = TRUE, nchains = NULL, 
gamma = 1, sdepsilon = 0.05, save = FALSE)
}
\arguments{
  \item{model}{
formula; <graph> ~ <model specifications>
}
\item{burn.in}{
count; number of burn-in iterations for each chain
}
  \item{main.iter}{
count; number of iterations for the main chain(s)
}
  \item{aux.iter}{
count; number of auxiliary iterations
}
  \item{sdprop}{
vector; standard deviation for the (independent) Normal proposals (for single-site update)
}
  \item{sdprior}{
vector; standard deviation for the (independent) Normal priors
}
  \item{mprior}{
vector; mean for the Normal priors
}
  \item{popMCMC}{
logical; If TRUE, population MCMC is performed using snooker update
}
  \item{nchains}{
count; number of chains for population MCMC
}
  \item{gamma}{
scalar; snooker factor
}
  \item{sdepsilon}{
parameter; standard deviation for the snooker update parameter (epsilon)
}
  \item{save}{
logical; if TRUE a file called "bergm.out" is saved in the working directory
}
  \item{...}{
additional arguments, to be passed to lower-level functions
}
}
\value{
\item{theta}{matrix containing all the parameter values in the chain(s)}
\item{dim}{model dimensions}
\item{chains}{number of chains}
\item{iter}{number of iterations for the main chain(s)}
\item{rate}{acceptance rate(s)}
\item{mod}{model specifications}
}
\seealso{
 \code{\link{mcmc.output}}, 
 \code{\link{bgof}}
}
\examples{
# load the Florentine marriage network 
# included with the ergm package

data(florentine)

# Estimation of a 2-dimensional model
# measuring the propensity to form 2-stars
# popMCMC with snooker update is used
# (this  will take a few minutes)

flo <- bergm(flomarriage~edges+kstar(2),
main.iter=5000,sdprior=c(5,10),gamma=0.8)

# MCMC diagnostics for the overall posterior
# and second chain

mcmc.output(flo,lags=200,chain=2)

# Bayesian goodness-of-fit test

bgof(flo,lags=150,n.sim=100,n.deg=10,n.dist=9,n.esp=6)
}