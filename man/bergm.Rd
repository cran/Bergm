\name{bergm}
\Rdversion{1.4}
\alias{bergm}
\title{
Bayesian exponential random graph models
}
\description{
Function to fit Bayesian exponential random graphs models
using the exchange algorithm with either single-site update
or population MCMC with parallel ADS move.
}
\usage{
bergm(model, 
      burn.in = 0, 
      main.iter = 5000, 
      aux.iter = 1000, 
      sdprop = NULL, 
      sdprior = NULL, 
      mprior = NULL, 
      popMCMC = TRUE, 
      nchains = NULL, 
      gamma = 1, 
      sdepsilon = 0.05, 
      save = FALSE,
      ...)
}
\arguments{
  \item{model}{
formula; an \code{R} formula object, of the form  <graph> ~ <model specifications> where <graph> is a \code{\link{network}} object and <model specifications> are \code{\link{ergm-terms}}.
}
\item{burn.in}{
count; number of burn-in iterations for the beginning of an MCMC run. If population MCMC is performed, it refers to the number of burn-in iterations for each chain of the population.
}
  \item{main.iter}{
count; number of iterations for the main chain(s).
}
  \item{aux.iter}{
count; number of auxiliary iterations used to simulate networks.
}
  \item{sdprop}{
vector; standard deviation of the (independent) Normal proposals (for single-site update procedure).
}
  \item{sdprior}{
vector; standard deviation of the (independent) Normal priors.
}
  \item{mprior}{
vector; mean of the Normal priors.
}
  \item{popMCMC}{
logical; If TRUE, population MCMC is performed using parallel ADS procedure.
}
  \item{nchains}{
count; number of chains for population MCMC procedure.
}
  \item{gamma}{
scalar; ``snooker move'' factor.
}
  \item{sdepsilon}{
parameter; standard deviation for ``snooker move'' parameter (epsilon).
}
  \item{save}{
logical; if TRUE a file called "bergm.out" is saved in the working directory.
}
  \item{...}{
additional arguments, to be passed to lower-level functions.
}
}
\value{
\item{theta}{matrix containing all the parameter values in the chain(s).}
\item{dim}{model dimensions.}
\item{chains}{number of chains.}
\item{iter}{number of iterations for the main chain(s).}
\item{rate}{acceptance rate(s).}
\item{mod}{model specifications.}
}
\seealso{
 \code{\link{mcmc.output}}, 
 \code{\link{posterior.plot}},
 \code{\link{bgof}}.
}
\examples{
# load the Florentine marriage network 
# included with the ergm package

data(florentine)

# Estimation of a 3-dimensional model
# measuring the propensity to form 2- and 3- stars.
# PopMCMC with parallel ADS approach is used
# (this  will take about 6 minutes)

flo <- bergm(flomarriage~edges+kstar(2:3),
burn.in=500,aux.iter=3000,main.iter=5000,
sdprior=c(5,5),gamma=0.6,sdepsilon=0.1)

# MCMC diagnostics for the second chain
# and overall poterior estimate

diagnostics <- mcmc.output(flo,lags=300,chain=2)

# Posterior density scatterplots 
# and covariance matrix

posterior.plot(diagnostics)

# Bayesian goodness-of-fit test

bgof(flo,lags=200,n.sim=100,n.deg=10,n.dist=9,n.esp=6)
}