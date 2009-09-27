\name{mcmc.output}
\Rdversion{1.1}
\alias{mcmc.output}
\title{
MCMC output diagnostics
}
\description{
This function returns the posterior density estimate and creates simple diagnostic plots for the MCMC produced from a fit. 
}
\usage{
mcmc.output(out, lags = 100, chain = 1, posterior = TRUE)
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
  \item{posterior}{
logical; if TRUE the MCMC output for the overall posterior density
}
}
\seealso{
  \code{\link{bergm}}
}
