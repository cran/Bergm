bergm <- function (formula, 
                   burn.in=100,
                   main.iters=1000,
                   aux.iters=1000, 
                   m.prior = NULL, 
                   sigma.prior = NULL, 
                   nchains = NULL, 
                   gamma = 0.5, 
                   sigma.epsilon = NULL,
                   save = FALSE,
                   ...){ 	
ptm=proc.time()   

y <- ergm.getnetwork(formula)
model <- ergm.getmodel(formula,y)  
Clist <- ergm.Cprepare(y, model)
Clist$maxedges = 1 + max(20000,Clist$nedges)
MHproposal <- MHproposal.ergm(formula,
                              constraints=~., 
                              arguments = NULL,
                              nw = y, 
                              model = model, 
                              weights = "default", #TNT sampler
                              class = "c")  
stats <- .C("network_stats_wrapper",
           as.integer(Clist$tails), as.integer(Clist$heads), 
           as.integer(Clist$nedges),
           as.integer(Clist$n),
           as.integer(Clist$dir), as.integer(Clist$bipartite), 
           as.integer(Clist$nterms), 
           as.character(Clist$fnamestring), as.character(Clist$snamestring), 
           as.double(Clist$inputs),
           gs = double(Clist$nstats),
           PACKAGE="ergm")$gs
snooker <- 0
if (is.null(m.prior))  m.prior <- rep(0,Clist$nstats) 
if (is.null(sigma.prior)) sigma.prior <- diag(100, Clist$nstats) 
if (is.null(nchains)) nchains <- 2*Clist$nstats
if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025,Clist$nstats)
if (Clist$nstats==1){ 
	nchains <- 1
	sigma.epsilon <- diag(gamma,Clist$nstats)	
}
Theta<- array(NA,c(main.iters,Clist$nstats,nchains))
theta <- matrix(runif(Clist$nstats*nchains,min=-.1,max=.1),
                Clist$nstats,nchains)
acc.counts <- rep(0, nchains)
pr <- theta1 <- rep(0,Clist$nstats)

tot.iters <- burn.in + main.iters

for (k in 1:tot.iters) {		
		for (h in 1:nchains) {
			if(Clist$nstats>1 && nchains>1){
				snooker <- gamma*apply(theta[,sample(seq(1,nchains)[-h],2)],1,diff)
			}
			
			theta1 <- theta[,h] + snooker + rmvnorm(1,sigma=sigma.epsilon)[1,]
			pr <- dmvnorm(rbind(theta1,theta[,h]),mean=m.prior,sigma=sigma.prior)
			prr <- pr[1]/pr[2]

			delta <- .C("MCMC_wrapper", as.integer(length(c(Clist$nedges, 0, 0))), 
			as.integer(c(Clist$nedges, 0, 0)), 
			as.integer(Clist$tails), as.integer(Clist$heads), 
			as.integer(Clist$maxpossibleedges), 
			as.integer(Clist$n), as.integer(Clist$dir), as.integer(Clist$bipartite), 
			as.integer(Clist$nterms), as.character(Clist$fnamestring), 
			as.character(Clist$snamestring), as.character(MHproposal$name), 
			as.character(MHproposal$package), as.double(Clist$inputs), 
			as.double(theta1), #
			as.integer(1), statsmatrix = double(1 * Clist$nstats), 
			as.integer(aux.iters), as.integer(0), 
			newnwtails = integer(Clist$maxedges), newnwheads = integer(Clist$maxedges), 
			as.integer(FALSE), as.integer(MHproposal$bd$attribs), 
			as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin), 
			as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin), 
			as.integer(MHproposal$bd$condAllDegExact), 
			as.integer(length(MHproposal$bd$attribs)), 
			as.integer(Clist$maxedges), PACKAGE="ergm")$statsmatrix 	

			beta <- t(theta[,h] - theta1) %*% delta + log(prr)
                
			if (beta >= log(runif(1))) {
				theta[,h] <- theta1
				if (k > burn.in) acc.counts[h] <- acc.counts[h] + 1
			}     
		}
		if (k > burn.in) Theta[k-burn.in,,] <- theta
}
if(nchains==1){
	Theta <- as.matrix(Theta[,,1])
}

time = proc.time() - ptm
    
out=list(Clist=Clist,
         MHproposal=MHproposal,
         formula=formula,
         model=model,
         nnodes=Clist$n,#
         specs=model$coef.names,
         dim=Clist$nstats,#
         nchains=nchains,
         stats=stats,
         Theta=Theta,
         nchains=nchains,
         acc.rates=acc.counts/main.iters,
         m.prior=m.prior,
         sigma.prior=sigma.prior,
         aux.iters=aux.iters,
         time=time) 
if (save == TRUE) dput(out, "bergm.out")
out
}

