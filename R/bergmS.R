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
	
	#print(post[[i]]$acc.rates)	
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
 
 
 
 #eta1 <- ergm.eta(theta1, m$etamap)       # non-CERGMs
z <- .C("MCMC_wrapper", as.integer(length(c(post[[m1]]$Clist$nedges, 0, 0))), 
as.integer(c(post[[m1]]$Clist$nedges, 0, 0)), 
as.integer(post[[m1]]$Clist$tails), as.integer(post[[m1]]$Clist$heads), as.integer(post[[m1]]$Clist$maxpossibleedges), 
as.integer(post[[m1]]$Clist$n), as.integer(post[[m1]]$Clist$dir), as.integer(post[[m1]]$Clist$bipartite), 
as.integer(post[[m1]]$Clist$nterms), as.character(post[[m1]]$Clist$fnamestring), 
as.character(post[[m1]]$Clist$snamestring), as.character(post[[m1]]$MHproposal$name), 
as.character(post[[m1]]$MHproposal$package), as.double(post[[m1]]$Clist$inputs), 
as.double(theta1), #
as.integer(1), 
statsmatrix = double(1 * post[[m1]]$Clist$nstats), 
as.integer(aux.iters), as.integer(0), 
newnwtails = integer(post[[m1]]$Clist$maxedges), newnwheads = integer(post[[m1]]$Clist$maxedges), 
as.integer(FALSE), as.integer(post[[m1]]$MHproposal$bd$attribs), 
as.integer(post[[m1]]$MHproposal$bd$maxout), as.integer(post[[m1]]$MHproposal$bd$maxin), 
as.integer(post[[m1]]$MHproposal$bd$minout), as.integer(post[[m1]]$MHproposal$bd$minin), 
as.integer(post[[m1]]$MHproposal$bd$condAllDegExact), as.integer(length(post[[m1]]$MHproposal$bd$attribs)), 
as.integer(post[[m1]]$Clist$maxedges), PACKAGE="ergm") 

nedges <- z$newnwtails[1]
if (nedges >= post[[m1]]$Clist$maxedges) {
	post[[m1]]$Clist$maxedges <- post[[m1]]$Clist$maxedges * 10
	z <- .C("MCMC_wrapper", as.integer(length(c(post[[m1]]$Clist$nedges, 0, 0))), 
	as.integer(c(post[[m1]]$Clist$nedges, 0, 0)), 
	as.integer(post[[m1]]$Clist$tails), as.integer(post[[m1]]$Clist$heads), as.integer(post[[m1]]$Clist$maxpossibleedges), 
	as.integer(post[[m1]]$Clist$n), as.integer(post[[m1]]$Clist$dir), as.integer(post[[m1]]$Clist$bipartite), 
	as.integer(post[[m1]]$Clist$nterms), as.character(post[[m1]]$Clist$fnamestring), 
	as.character(post[[m1]]$Clist$snamestring), as.character(post[[m1]]$MHproposal$name), 
	as.character(post[[m1]]$MHproposal$package), as.double(post[[m1]]$Clist$inputs), 
	as.double(theta1), #
	as.integer(1), 
	statsmatrix = double(1 * post[[m1]]$Clist$nstats), 
	as.integer(aux.iters), as.integer(0), 
	newnwtails = integer(post[[m1]]$Clist$maxedges), newnwheads = integer(post[[m1]]$Clist$maxedges), 
	as.integer(FALSE), as.integer(post[[m1]]$MHproposal$bd$attribs), 
	as.integer(post[[m1]]$MHproposal$bd$maxout), as.integer(post[[m1]]$MHproposal$bd$maxin), 
	as.integer(post[[m1]]$MHproposal$bd$minout), as.integer(post[[m1]]$MHproposal$bd$minin), 
	as.integer(post[[m1]]$MHproposal$bd$condAllDegExact), as.integer(length(post[[m1]]$MHproposal$bd$attribs)), 
	as.integer(post[[m1]]$Clist$maxedges), PACKAGE="ergm") 
}else if (nedges == 0) {
	newedgelist <- matrix(0, ncol = 2, nrow = 0)
}else {
	newedgelist <- cbind(z$newnwtails[2:(nedges + 1)], z$newnwheads[2:(nedges + 1)])
}
nw <- network.update(y, newedgelist, matrix.type = "edgelist")
nw <- as.network.uncompressed(nw)

# update Clist for model m2		
Clist <- ergm.Cprepare(nw, post[[m]]$model)
Clist$maxedges = 1 + max(20000,Clist$nedges)	
mstats <- .C("network_stats_wrapper",
             as.integer(Clist$tails), as.integer(Clist$heads), 
             as.integer(Clist$nedges),as.integer(Clist$n),
             as.integer(Clist$dir), as.integer(Clist$bipartite), 
             as.integer(Clist$nterms), as.character(Clist$fnamestring), 
             as.character(Clist$snamestring), as.double(Clist$inputs),
             gs = double(Clist$nstats),PACKAGE="ergm")$gs

deltad <- z$statsmatrix
     
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

