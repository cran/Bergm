\name{mcmc.output}
\Rdversion{1.4}
\alias{mcmc.output}
\title{
MCMC output diagnostics
}
\description{
This function returns the posterior density estimate and creates simple diagnostic plots for the MCMC produced from a fit. 
}
\usage{
mcmc.output(out, 
            lags = 100, 
            chain = 0)
}
\arguments{
  \item{out}{ 
  an \code{R} object of class \code{bergm}.
}
  \item{lags}{
count; number of lags for the autocorrelation plots.
}
  \item{chain}{
count; used to plot the MCMC output of the selected chain (only for population MCMC).
}
}
\value{
\item{posterior}{matrix containing the MCMC output for the main chain.}
\item{dim}{model dimensions.}
}
\seealso{
  \code{\link{bergm}}.
}
