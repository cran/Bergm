% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evidence.R
\name{evidence}
\alias{evidence}
\title{Wrapper function for evidence estimation}
\usage{
evidence(evidence.method = c("CJ", "PP"), ...)
}
\arguments{
\item{evidence.method}{vector Method to estimate the marginal likelihood. Options are: \code{"CJ"}, in
which case the marginal likelihood is estimated with Chib and Jeliazkov's method; \code{"PP"}, in
which case the marginal likelihood is estimated with Power posteriors.}

\item{...}{further arguments to be passed.
See \code{evidenceCJ} and \code{evidencePP}.}
}
\description{
Function to estimate the evidence (marginal likelihood) with Chib and Jeliazkov's method
or Power posteriors, based on the adjusted pseudolikelihood function.
}
\examples{
\dontrun{
# Load the florentine marriage network:
data(florentine)
                                                
# MCMC sampling and evidence estimation:
CJE <- evidence(evidence.method = "CJ",
                formula     = flomarriage ~ edges + kstar(2),
                main.iters  = 30000,
                burn.in     = 2000,
                aux.iters   = 1000,
                num.samples = 25000,
                V.proposal  = 2.5,
                ladder      = 100,
                seed        = 1)
                                   
# Posterior summaries:
summary(CJE)

# MCMC diagnostics plots:
plot(CJE)
    
# Log-evidence (marginal likelihood) estimate:
CJE$log.evidence
}

}
\references{
Bouranis, L., Friel, N., & Maire, F. (2018). Bayesian model selection for exponential 
random graph models via adjusted pseudolikelihoods. 
Journal of Computational and Graphical Statistics, 27(3), 516-528. 
\url{https://arxiv.org/abs/1706.06344}
}
