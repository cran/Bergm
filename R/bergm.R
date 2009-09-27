bergm <-
function(model,
		 burn.in=1000,main.iter=25000,
		 sdprop=NULL,sdprior=50,mprior=0,theta=NULL,
		 popMCMC=FALSE,
		 nchains=NULL,block.iter=1000,sdblock=NULL,
		 sdgamma=0.5,sdepsilon=0.05,
		 save=FALSE){
				  	          
	mod <- ergm.getmodel(model,ergm.getnetwork(model))
	stat <- ergm.getglobalstats(ergm.getnetwork(model),mod)
	mrow <- length(stat)
	
	if(popMCMC==FALSE){ 
		mcol <- mrow
		block.iter=0
		H <- matrix(0,main.iter,mrow)
		if(is.null(theta)) theta <- rep(0,mrow)
		if(is.null(sdprop)) sdprop <- rep(0.1,mrow)
            nchains <- 1
	}else{ 
		if(is.null(nchains)) nchains <- mrow+1
		mcol <- nchains
		H <- matrix(0,main.iter,mrow*mcol)
		if(is.null(theta)){theta <- matrix(0,mrow,mcol)
		}else{theta <- matrix(theta,mrow,mcol,byrow=FALSE)}
		if(is.null(sdblock)) sdblock <- rep(0.1,mrow)
	}
	
	iter <- block.iter+main.iter
	pr <- rep(0,mrow)
 	thetad <- theta
	accept <- rep(0,mcol)

	for(k in 1:iter){
		for(h in 1:mcol){
			
			# Single-site update
			if(popMCMC==FALSE){
		        	thetad[h] <- rnorm(1,theta[h],sdprop[h])
				pr <- dnorm(c(theta,thetad),mprior,sdprior)
				prr <- prod(pr[1:mrow]/pr[1:mrow])    
				yd <- simulate(model,
       	         		       theta0 = thetad, 
                       		       burnin = burn.in)
				delta <- ergm.getglobalstats(yd,mod)-stat
				beta <- t(theta-thetad)%*%delta+log(prr)
				if(beta >= log(runif(1))){ 
					theta[h] <- thetad[h]; accept[h] <- accept[h]+1
				}
			}
			
			# PopMCMC	
			else{
			# Block update
    			if(k <= block.iter){ 
    				thetad[,h] <- rnorm(mrow,theta[,h],sdblock)
				}
			# Snooker update
    			else{
    				if(h == mcol) thetad[,h] <- theta[,h]+rnorm(1,0,sdgamma)*
    				     	      (theta[,h]-theta[,1])+rnorm(mrow,0,sdepsilon)     
				else thetad[,h] <- theta[,h]+rnorm(1,0,sdgamma)*
                                                   (theta[,h]-theta[,h+1])+rnorm(mrow,0,sdepsilon)
    			}	       
			pr <- dnorm(c(theta[,h],thetad[,h]),mprior,sdprior)
			prr <- prod(pr[1:mrow]/pr[1:mrow])   
			yd <- simulate(model,
                                       theta0 = thetad[,h], 
                                       burnin = burn.in)
			delta <- ergm.getglobalstats(yd,mod)-stat
			beta <- t(theta[,h]-thetad[,h])%*%delta+log(prr)
				if(beta >= log(runif(1))){ 
					theta[,h] <- thetad[,h] 
					if(k > block.iter){accept[h] <- accept[h]+1}
				}	
			}	
		}	
		if(k > block.iter) H[k-block.iter,] <- theta
	}
        out = list(theta=H,dim=mrow,chains=nchains,iter=main.iter,rate=(accept/main.iter),mod=model)
        if(save==TRUE) dput(out,"bergm.out")
        out
}

