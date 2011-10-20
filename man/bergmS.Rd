\name{bergmS}
\alias{bergmS}
\title{
Bayesian model selection for exponential random graph models
}
\description{
Bayesian model selection for exponential random graphs models using the auto-RJ exchange algorithm.
}
\usage{
bergmS(formulae, 
       iters = 10000, 
       m.priors = NULL, 
       sigma.priors = NULL, 
       gammas = NULL, 
       nchains = NULL, 
       sigma.epsilons = NULL, 
       aux.iters = 1000, 
       main.iters = NULL,
       burn.ins=NULL,
       save=FALSE, 
       ...)
}
\arguments{
  \item{formulae}{
list; \code{R} formula objects (competing models) of the form  <graph> ~ <model specifications>, see \code{\link{ergm-terms}}.
}
  \item{iters}{
count; number of iterations for auto-RJ exchange algorithm (online step).
}
  \item{m.priors}{
list; means of the multivariate Normal priors of each competing model.
}
  \item{sigma.priors}{
list; means of the multivariate Normal priors of each competing model.
}
  \item{gammas}{
vector; ``snooker move factors'' for each competing model (offline step).
}
  \item{nchains}{
vector; number of chains for each competing model (offline step).
}
  \item{sigma.epsilons}{
list; containing variance/covariance matrices for the multivariate Normal proposals or ``snooker move parameters'' for each competing model (offline step).
}
  \item{aux.iters}{
count; number of auxiliary iterations for network simulation for each competing model.
}
  \item{main.iters}{
vector; number of iterations for the main chain(s) for each competing model (offline step).
}
  \item{burn.ins}{
vector; number of burn-in iterations for the beginning of an MCMC run for each competing model (offline step).
}
  \item{save}{
logical; if TRUE a file called "bergmS.out" is saved in the working directory.
}
  \item{\dots}{
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
 \code{\link{bergm}},
}
\examples{
# load the Florentine marriage network 
# included with the ergm package

data(florentine)
y <- flomarriage

# Competing models:

formulae <- c(y ~ edges, 
              y ~ edges + gwesp(log(2),fixed=TRUE),
              y ~ edges + gwdegree(log(2),fixed=TRUE))

# Model selection via auto-RJ exchange algorithm
# (this  will take about 3.5 minutes)

flo <- bergmS(formulae,
                  iters=10000,
                  aux.iters=3000,
                  main.iters=rep(1000,3),
                  burn.in=rep(100,3),
                  gammas=c(0.8,1,1))

# MCMC diagnostics

flo.out <- bergmS.output(flo)

# Bayesian goodness-of-fit test

bgof(flo.out,lags=200,n.sim=100,n.deg=10,n.dist=10,n.esp=5)

}