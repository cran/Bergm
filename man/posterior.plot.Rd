\name{posterior.plot}
\Rdversion{1.4}
\alias{posterior.plot}
\title{
Posterior density scatterplots
}
\description{
This function produces two-dimensional scatterplots of the estimated posterior density and returns the variance/covariance matrix of the parameters.
}
\usage{
posterior.plot(mcmc.out,
               ...)
}
\arguments{
  \item{mcmc.out}{ 
  an \code{R} object of class \code{mcmc.output}.
}
  \item{...}{
additional arguments, to be passed to lower-level functions.
}
}
\seealso{
  \code{\link{bergm}}.
}
