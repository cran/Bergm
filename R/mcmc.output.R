mcmc.output <- function(out,lags=100,chain=1,save=FALSE){

cat('Post. Mean:',"\n")
postmean <- matrix(apply(as.matrix(out$theta[,]),2,mean),out$dim,out$chains,byrow=FALSE)
colnames(postmean) <- paste("chain",seq(1,out$chains))
rownames(postmean) <- paste("theta",seq(1,out$dim))
postmean <- as.table(postmean)
print(postmean)

cat("\n",'Post. Sd.:',"\n")
postsd <- matrix(apply(as.matrix(out$theta[,]),2,sd),out$dim,out$chains,byrow=FALSE)
colnames(postsd) <- paste("chain",seq(1,out$chains))
rownames(postsd) <- paste("theta",seq(1,out$dim))
postsd <- as.table(postsd)
print(postsd)

cat("\n")
rates <- matrix(out$rate,1,out$chains)
colnames(rates) <- paste("chain",seq(1,out$chains))
rownames(rates) <- paste("Acc. rates:")
rates <- as.table(rates)
print(rates)

F <- as.matrix(out$theta[,1:out$dim])
if(out$chains >= 2){
	for(h in 2:out$chains){
		F <- rbind(F,out$theta[,((out$dim*(h-1)+1)):(out$dim*h)])
	}
}
cat("\n",'Overall Posterior Density Estimates:',"\n")
overall <- matrix(c(apply(F,2,mean),apply(F,2,sd)),out$dim,2)
colnames(overall) <- c("Post. Mean:","Post. Sd.:")
rownames(overall) <- paste("theta",seq(1,out$dim))
all <- as.table(overall)
print(overall)

par(mfrow=c(1,out$dim),oma=c(0,0,4,0))
for(i in 1:out$dim){
plot(density(F[,i]),
main="",
axes=FALSE, xlab=list(bquote(theta[.(i)]),cex=1.5),ylab="",lwd=2)
axis(1); axis(2)
}
title("Posterior density",outer=TRUE)

dev.new()

    G<-mcmc(data=out$theta,start=1,thin=1)
	par(mfrow=c(out$dim,3),oma=c(0,0,4,0))
	for(i in (((chain-1)*out$dim)+1):(out$dim*chain)){
		if(chain >= 2){j <- i-(out$dim*(chain-1))}else{j <- i}
		
		plot(density(out$theta[,i]),
		main="",
		axes=FALSE, xlab=list(bquote(theta[.(j)]),cex=1.5),ylab="",lwd=2)
		axis(1); axis(2)

		plot(out$theta[,i],type="l",xlab="iterations",ylab="")
		autocorr.plot(G[,i],lags,auto.layout=FALSE)
	}	
	title(paste("MCMC output ( chain =",chain,")"),outer=TRUE)
    out = list(postmean=postmean,postsd=postsd,accrate=rates,all=overall)
    if(save==TRUE) dput(out,"mcmc.output.out")
}