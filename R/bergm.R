bergm <- function (model, burn.in = 0, main.iter = 5000, aux.iter = 1000, 
    sdprop = NULL, sdprior = NULL, mprior = NULL, 
    popMCMC = TRUE, nchains = NULL, gamma = 1, sdepsilon = 0.05, 
    save = FALSE) 
{
    mod <- ergm.getmodel(model, ergm.getnetwork(model))
    stat <- ergm.getglobalstats(ergm.getnetwork(model), mod)
    mrow <- length(stat)
    if (is.null(sdprior)) {
        sdprior <- rep(10, mrow)
    }
    if (is.null(mprior)) {
        mprior <- rep(0, mrow)
    }
    if (popMCMC == FALSE) {
        mcol <- mrow
        H <- matrix(0, main.iter, mrow)
        theta <- runif(mrow, max = 0.1)
        if (is.null(sdprop)) {
            sdprop <- rep(0.1, mrow)
        }
        nchains <- 1
    }
    else {
        if (is.null(nchains)) {
            nchains <- 2 * mrow
        }
        mcol <- nchains
        H <- matrix(0, main.iter, mrow * mcol)
        theta <- matrix(runif(mrow * mcol, max = 0.1), mrow, mcol)
    }
    iter <- burn.in + main.iter
    pr <- rep(0, mrow)
    thetad <- rep(0, mrow)
    accept <- rep(0, mcol)
    for (k in 1:iter) {
        for (h in 1:mcol) {
            if (popMCMC == FALSE) {
                thetad[h] <- rnorm(1, theta[h], sdprop[h])
                pr <- dnorm(theta, mprior, sdprior)
                prd <- dnorm(thetad, mprior, sdprior)
                prr <- prod(prd/pr)
                yd <- simulate(model, theta0 = thetad, burnin = aux.iter)
                delta <- ergm.getglobalstats(yd, mod) - stat
                beta <- t(theta - thetad) %*% delta + log(prr)
                if (beta >= log(runif(1))) {
                  theta[h] <- thetad[h]
                  if (k > burn.in) {
                    accept[h] <- accept[h] + 1
                  }
                }
            }
            else {
                thetad <- theta[, h] + 
                   gamma * apply(theta[, sample(seq(1, mcol)[-h], 2)], 1, diff) + 
                   rnorm(mrow, 0, sdepsilon)
                pr <- dnorm(theta[, h], mprior, sdprior)
                prd <- dnorm(thetad, mprior, sdprior)
                prr <- prod(prd/pr)
                yd <- simulate(model, theta0 = thetad, burnin = aux.iter)
                delta <- ergm.getglobalstats(yd, mod) - stat
                beta <- t(theta[, h] - thetad) %*% delta + log(prr)
                if (beta >= log(runif(1))) {
                  theta[, h] <- thetad
                  if (k > burn.in) {
                    accept[h] <- accept[h] + 1
                  }
                }
            }
        }
        if (k > burn.in) 
            H[k - burn.in, ] <- theta
    }
    out = list(theta = H, dim = mrow, chains = nchains, iter = main.iter, 
        rate = (accept/main.iter), mod = model)
    if (save == TRUE) {
        dput(out, "bergm.out")
    }
    out
}
