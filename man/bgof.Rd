\name{bgof}
\Rdversion{1.1}
\alias{bgof}
\title{
Bayesian goodness-of-fit diagnostics
}
\description{
Calculates summaries for degree, minimum geodesic distances, and edge-wise shared partner distributions to diagnose the Bayesian goodness-of-fit of exponential family random graph models.
}
\usage{
bgof(out, directed = FALSE, lag = 100, n.sim = NULL, burn.in = 10000, n.deg = NULL, n.dist = NULL, n.esp = NULL, n.ideg = NULL, n.odeg = NULL, save = FALSE)
}
\arguments{
  \item{out}{
\code{bergm} output
}
  \item{directed}{
logical; TRUE if the observed graph is directed
}
  \item{lag}{
count; interval between each sampled posterior parameter estimate
}
  \item{n.sim}{
count; number of sampled parameter values 
}
  \item{burn.in}{
count; number of iterations used for each graph simulation
}
  \item{n.deg}{
count; used to plot only the first \code{n.deg-1} degree distributions
}
  \item{n.dist}{
count; used to plot only the first \code{n.dist-1} geodesic distances distributions
}
  \item{n.esp}{
count; used to plot only the first \code{n.esp-1} edge-wise shared partner distributions
}
  \item{n.ideg}{
count; used to plot only the first \code{n.ideg-1} in-degree distributions
}
  \item{n.odeg}{
count; used to plot only the first \code{n.odeg-1} out-degree distributions
}
  \item{save}{
logical; if TRUE a file called "bgof.out" is saved in the working directory
}
}
\examples{
data(florentine)

flo1 <- bergm(flobusiness~edges+kstar(2),main.iter=15000,
sdprop=c(1,0.1),sdprior=20,theta=c(-2,0.1))

## Bayesian goodness-of-fit test 

bgof(flo1,lag=150,n.sim=100)

## Let's improve the plots

bgof(flo1,lag=150,n.sim=100,n.deg=10,n.dist=9,n.esp=6)

}
