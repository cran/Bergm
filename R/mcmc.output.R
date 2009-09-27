mcmc.output <- function(out,lags=100,chain=1,posterior=TRUE){

	if(out$chains > 1){
	
	cat("\n",'Posterior mean:',"\n")
	postmean <- matrix(apply(as.matrix(out$theta[,]),2,mean),
					   out$dim,out$chains,byrow=FALSE)
	colnames(postmean) <- paste("Chain",seq(1,out$chains))
	rownames(postmean) <- paste("theta",seq(1,out$dim))
	postmean <- as.table(postmean)
	print(postmean)

	cat("\n",'Posterior standard deviation:',"\n")
	postsd <- matrix(apply(as.matrix(out$theta[,]),2,sd),
					 out$dim,out$chains,byrow=FALSE)
	colnames(postsd) <- paste("Chain",seq(1,out$chains))
	rownames(postsd) <- paste("theta",seq(1,out$dim))
	postsd <- as.table(postsd)
	print(postsd)

	cat("\n")
	rates <- matrix(out$rate,1,out$chains)
	colnames(rates) <- paste("Chain",seq(1,out$chains))
	rownames(rates) <- paste("Acceptance rates:")
	rates <- as.table(rates)
	print(rates)
	
		F <- as.matrix(out$theta[,1:out$dim])
			for(h in 2:out$chains){
				F <- rbind(F,out$theta[,((out$dim*(h-1)+1)):(out$dim*h)])
			}
		cat("\n",'Overall Posterior Density Estimates:',"\n")
		overall <- matrix(c(apply(F,2,mean),apply(F,2,sd)),out$dim,2)
		colnames(overall) <- c("Post. mean:","Post. sd.:")
		rownames(overall) <- paste("theta",seq(1,out$dim))
		all <- as.table(overall)
		print(overall)

		if(posterior==TRUE){
			K <- mcmc(data=F,start=1,thin=1)
			par(mfrow=c(out$dim,3),oma=c(0,0,4,0))
				for(i in 1:out$dim){
					plot(density(F[,i]),
						 main="",
						 axes=FALSE, 
						 xlab=list(bquote(theta[.(i)]),cex=1.5),
						 ylab="",lwd=2)
					axis(1); axis(2)
					plot(F[,i],type="l",xlab="Iterations",ylab="")
					autocorr.plot(K[,i],lags,auto.layout=FALSE)
				}
			title("Overall posterior estimate",outer=TRUE)
		}
	}else{
		
		F <- out$theta
		cat("\n",'Posterior Density Estimates:',"\n")
		overall <- matrix(c(apply(F,2,mean),apply(F,2,sd)),out$dim,2)
		colnames(overall) <- c("Post. mean:","Post. sd.:")
		rownames(overall) <- paste("theta",seq(1,out$dim))
		all <- as.table(overall)
		print(overall)	
		
		cat("\n")
	  	rates <- matrix(out$rate,1,out$chains)
		colnames(rates) <- paste("")
		rownames(rates) <- paste("Acceptance rate:")
		rates <- as.table(rates)
		print(rates)	
		
	}
	

dev.new()

	
    G <- mcmc(data=out$theta,start=1,thin=1)
	par(mfrow=c(out$dim,3),oma=c(0,0,4,0))
	for(i in (((chain-1)*out$dim)+1):(out$dim*chain)){
		if(chain >= 2){ j <- i-(out$dim*(chain-1)) }else{ j <- i }
		plot(density(out$theta[,i]),
			 main="",
			 axes=FALSE, 
			 xlab=list(bquote(theta[.(j)]),cex=1.5),ylab="",lwd=2)
		axis(1); axis(2)
		plot(out$theta[,i],type="l",xlab="Iterations",ylab="")
		autocorr.plot(G[,i],lags,auto.layout=FALSE)
	}	
	title(paste("MCMC output ( chain =",chain,")"),outer=TRUE)
}