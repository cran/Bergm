\name{bergm}
\Rdversion{2.0}
\alias{bergm}
\title{
Bayesian exponential random graph models
}
\description{
Function to fit Bayesian exponential random graphs models using the exchange algorithm with either single-site update or population MCMC with parallel ADS move.
}
\usage{
bergm(formula, 
      burn.in=NULL,
      main.iters=NULL,
      aux.iters=1000, 
      m.prior = NULL, 
      sigma.prior = NULL, 
      nchains = NULL, 
      gamma = NULL, 
      sigma.epsilon = NULL,
      save = FALSE,
      ...)	
}
\arguments{
  \item{formula}{
formula; an \code{R} formula object, of the form  <graph> ~ <model specifications> where <graph> is a \code{\link{network}} object and <model specifications> are \code{\link{ergm-terms}}.
}
\item{burn.in}{
count; number of burn-in iterations for the beginning of an MCMC run. If population MCMC is performed, it refers to the number of burn-in iterations for each chain of the population.
}
  \item{main.iters}{
count; number of iterations for the main chain(s).
}
  \item{aux.iters}{
count; number of auxiliary iterations used to simulate networks.
}
  \item{m.prior}{
vector; mean of the multivariate Normal prior.
}
  \item{sigma.prior}{
variance/covariance matrix for the multivariate Normal prior.
}
  \item{nchains}{
count; number of chains for ADS procedure.
}
  \item{gamma}{
scalar; ``snooker move factor.'' In case of one-dimensional models, the ADS procedure is disabled and \code{gamma} will be used as variance of the Normal proposal distribution.
}
  \item{sigma.epsilon}{
variance/covariance matrix for the multivariate Normal proposal or ``snooker move parameter''. In case of one-dimensional models, if \code{gamma} > 0, \code{sigma.espilon} is set equal to \code{gamma}.
}
  \item{save}{
logical; if TRUE a file called "bergm.out" is saved in the working directory.
}
  \item{...}{
additional arguments, to be passed to lower-level functions.
}
}
\value{
\item{Theta}{matrix containing all the parameter values in the chain(s).}
\item{dim}{model dimensions.}
\item{nchains}{number of chains.}
\item{time}{CPU time and real time (in seconds) that elapsed since the execution of the function commenced.}
\item{acc.rates}{acceptance rate(s).}
\item{formula}{model specifications.}
\item{nnodes}{number of nodes of the observed network.}
}
\seealso{
 \code{\link{bergm.output}}, 
 \code{\link{bgof}}.
}
\examples{
# load the Florentine marriage network 
# included with the ergm package

data(florentine)

# Estimation of a 3-dimensional model
# measuring the propensity to form 2- and 3- stars.
# Parallel ADS approach is used
# (this  will take about 1 minute)

flo <- bergm(flomarriage ~ edges + kstar(2:3),
             burn.in=500,
             aux.iters=3000,
             main.iters=1500)

# MCMC diagnostics

bergm.output(flo)

# Bayesian goodness-of-fit test

bgof(flo,n.sim=100,n.deg=10,n.dist=9,n.esp=6)
}