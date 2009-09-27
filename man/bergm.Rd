\name{bergm}
\Rdversion{1.1}
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
bergm(model, burn.in = 1000, main.iter = 25000, sdprop = NULL, sdprior = 50, mprior = 0, theta = NULL, popMCMC = FALSE, nchains = NULL, block.iter = 1000, sdblock = NULL, sdgamma = 0.5, sdepsilon = 0.05, save = FALSE)
}
\arguments{
  \item{model}{
formula; <graph> ~ <model specifications>
}
\item{burn.in}{
count; number of auxiliary iterations
}
  \item{main.iter}{
count; number of iterations for the main chain(s)
}
  \item{sdprop}{
vector; standard deviation for the (Normal) proposal distribution
}
  \item{sdprior}{
vector; standard deviation for the (Normal) prior distribution
}
  \item{mprior}{
vector; mean for the (Normal) prior distribution
}
  \item{theta}{
vector; initial parameter value
}
  \item{popMCMC}{
logical; If TRUE, population MCMC is performed using snooker update
}
  \item{nchains}{
count; number of chains for population MCMC
}
  \item{block.iter}{
count; number of iterations for the block update used at the start
of population MCMC procedure
}
  \item{sdblock}{
vector; standard deviation for the block update proposal
}
  \item{sdgamma}{
vector; standard deviation for the first snooker update parameter (gamma)
}
  \item{sdepsilon}{
vector; standard deviation for the second snooker update parameter (epsilon)
}
  \item{save}{
logical; if TRUE a file called "bergm.out" is saved in the working directory
}
}
\value{
\item{theta}{matrix containing all the parameter values estimated}
\item{dim}{model dimensions}
\item{chains}{number of chains}
\item{iter}{number of iterations for the main chain(s)}
\item{rate}{acceptance rate(s)}
\item{mod}{model}
}

\references{
Caimo, A. and Friel, N. (2009) ``Bayesian inference for exponential random graph models,'' (submitted).

Murray, I., Ghahramani, Z., and MacKay, D. (2006), 
``MCMC for doubly-intractable distributions,'' in Proceedings of the 22nd Annual Conference on Uncertainty in Artificial Intelligence (UAI-06), Arlington, Virginia: AUAI Press.

Handcock, M. S., Hunter, D. R., Butts, C. T., Goodreau, S. M., and Morris, M. (2007), ``statnet: Software tools for the representation, visualization,
analysis and simulation of network data,'' Journal of Statistical Software, 24,1-11.

Hunter, D. R., Handcock, M. S., Butts, C. T., Goodreau, S. M., and Morris M. (2008), ``ergm: A Package to Fit, Simulate and Diagnose Exponential-Family Models for Networks'', 24, 1-29.
}

\examples{
data(florentine)

# Bayesian estimation via exchange algorithm
# (population MCMC with snooker update)

flo2 <- bergm(flobusiness~edges+kstar(2),main.iter=10000,           popMCMC=TRUE,sdprior=50,block.iter=1000,sdgamma=1)
}