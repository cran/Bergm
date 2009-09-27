\name{mcmc.output}
\Rdversion{1.1}
\alias{mcmc.output}
\title{
MCMC output diagnostics
}
\description{
This function returns the posterior density estimates and creates simple diagnostic plots for the MCMC produced from a fit. 
}
\usage{
mcmc.output(out, lags = 100, chain = 1, save = FALSE)
}
\arguments{
  \item{out}{ 
\code{bergm} output
}
  \item{lags}{
count; number of lags for the autocorrelation plots
}
  \item{chain}{
count; used to plot the MCMC output for the selected chain
}
  \item{save}{
logical; if TRUE a file called "mcmc.output.out" is saved in the working directory
}
}
\examples{
data(florentine)

flo2 <- bergm(flobusiness~edges+kstar(2),main.iter=8000,           popMCMC=TRUE,sdprior=50,block.iter=1000,sdgamma=1)

## MCMC output diagnostics for the second chain

mcmc.output(flo2,lags=150,chain=2)

}
