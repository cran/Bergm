bergmS <- function(formulae,
                         iters=10000,
                         m.priors=NULL,
                         sigma.priors=NULL,
                         gammas=NULL,
                         nchains=NULL,
                         sigma.epsilons=NULL,
                         aux.iters=1000,
                         main.iters=NULL,
                         burn.ins=NULL,
                         save=FALSE,
                         ...){

ptm=proc.time()             	
nmodels <- length(formulae)
dims <- rep(0,nmodels)
y <- ergm.getnetwork(formulae[[1]])

if(is.null(main.iters)) main.iters <- rep(1000,nmodels) 
if(is.null(burn.ins)) burn.ins <- rep(100,nmodels) 
if(is.null(gammas)) gammas <- rep(0.5,nmodels)   	 
models <- specs <- theta <- Sigma <- Mu <- post <- stats <- vector("list", nmodels) 

for(i in 1:nmodels){ 	
	post[[i]] <- bergm(formula=formulae[[i]],
	                   burn.in=burn.ins[i],
	                   main.iters=main.iters[i],
	                   aux.iters=aux.iters,
	                   m.prior=m.priors[[i]],
	                   sigma.prior=sigma.priors[[i]],
	                   gamma=gammas[i],
	                   nchains=nchains[i],
	                   sigma.epsilon=sigma.epsilons[[i]])
	                   
	if(post[[i]]$nchains > 1){
		Sigma[[i]] <- cov(apply(post[[i]]$Theta,2,cbind))
		diag(Sigma[[i]])=8/7*diag(Sigma[[i]])
	}else{
		Sigma[[i]] <- var(post[[i]]$Theta)
		diag(Sigma[[i]])=8/7*Sigma[[i]]
	}
	Mu[[i]] <- apply(post[[i]]$Theta,2,mean)
	theta[[i]] <- runif(post[[i]]$dim,min=-0.1,max=0.1)
	dims[i] <- post[[i]]$dim	
	
	specs[[i]] <- post[[i]]$specs
}

Theta <- matrix(0,iters,max(dims))       
M <- rep(0,iters)
m <- sample(seq(1:nmodels),1)
Wrate <- matrix(0,2,nmodels)
Brate <- matrix(0,2,1)

for(i in 1:iters){

    m1 <- sample(1:nmodels,1) 
	if(m1!=m){
		Brate[2,] <- Brate[2,]+1
	}else{
		Wrate[2,m] <- Wrate[2,m]+1
	}
	theta1 <- rmvnorm(1,mean=Mu[[m1]],sigma=Sigma[[m1]])[1,]
    ww1 <- dmvnorm(theta1,mean=Mu[[m1]],sigma=Sigma[[m1]])
    pr1 <- dmvnorm(theta1,mean=post[[m1]]$m.prior,sigma=post[[m1]]$sigma.prior)    
    ww <- dmvnorm(theta[[m]],mean=Mu[[m]],sigma=Sigma[[m]])
    pr <- dmvnorm(theta[[m]],mean=post[[m]]$m.prior,sigma=post[[m]]$sigma.prior)
	prr <- pr1/pr
	wwr <- ww/ww1
 
nw <- simulate(formulae[[m1]], coef=theta1,
                  control=control.simulate.formula(MCMC.burnin=aux.iters,
                                                   MCMC.interval=0))
mstats <- summary(ergm.update.formula(formulae[[m]], nw~.))
deltad <- summary(ergm.update.formula(formulae[[m1]], nw~.)) - summary(formulae[[m1]])

    delta <- mstats - post[[m]]$stats
	beta <- t(theta[[m]]) %*% delta - t(theta1) %*% deltad + log(prr) + log(wwr)

	if (beta >= log(runif(1))){
		theta[[m1]] <- theta1
		if(m1!=m){
			Brate[1,] <- Brate[1,]+1
			m <- m1
		}else{
			Wrate[1,m] <- Wrate[1,m]+1
		}
	}
	Theta[i,1:dims[m]] <- theta[[m]]
	M[i] <- m
} 
time = proc.time()-ptm              

out=list(M=M,
         iters=iters,
         Theta=Theta,
         formulae=formulae,
         models=models,
         specs=specs,
         dims=dims,
         nmodels=nmodels,
         Baccept=Brate[1,]/Brate[2,],
         Waccept=Wrate[1,]/Wrate[2,],
         time=time)
if (save == TRUE) dput(out, "bergmS.out")
out
}

