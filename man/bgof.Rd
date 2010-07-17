\name{bgof}
\Rdversion{1.3}
\alias{bgof}
\title{
Bayesian goodness-of-fit diagnostics
}
\description{
Calculates summaries for degree, minimum geodesic distances, and edge-wise shared partner distributions to diagnose the Bayesian goodness-of-fit of exponential family random graph models.
}
\usage{
bgof(out, 
     directed = FALSE, 
     lags = 100, 
     n.sim = NULL, 
     aux.iter = 10000, 
     n.deg = NULL, 
     n.dist = NULL, 
     n.esp = NULL, 
     n.ideg = NULL, 
     n.odeg = NULL, 
     save = FALSE)
}
\arguments{
  \item{out}{
an \code{R} object of class \code{bergm}.
}
  \item{directed}{
logical; TRUE if the observed graph is directed.
}
  \item{lags}{
count; interval between each sampled posterior parameter value.
}
  \item{n.sim}{
count; number of sampled parameter values. 
}
  \item{aux.iter}{
count; number of iterations used for network simulation.
}
  \item{n.deg}{
count; used to plot only the first \code{n.deg-1} degree distributions.
}
  \item{n.dist}{
count; used to plot only the first \code{n.dist-1} geodesic distances distributions.
}
  \item{n.esp}{
count; used to plot only the first \code{n.esp-1} edge-wise shared partner distributions.
}
  \item{n.ideg}{
count; used to plot only the first \code{n.ideg-1} in-degree distributions.
}
  \item{n.odeg}{
count; used to plot only the first \code{n.odeg-1} out-degree distributions.
}
  \item{save}{
logical; if TRUE a file called "bgof.out" is saved in the working directory.
}
}
\seealso{
\code{\link{bergm}}.
}
