posterior.plot <- function(mcmc.out,...){
par(mfrow=c(mcmc.out$dim,mcmc.out$dim),mar=rep(1,4))
  for(i in 1:mcmc.out$dim){
    for(j in 1:i) plot(-1:1,-1:1,type="n",axes=FALSE,xlab="",ylab="")
    text(0,0,bquote(theta[.(i)]),cex=3)
      if(j < mcmc.out$dim){
        for(k in (i+1):mcmc.out$dim){
          smoothScatter(mcmc.out$posterior[,c(k,i)],
          nrpoints=0,nbin=200,...)
        }
      }
  }
  cat("\n",'Posterior covariance matrix:',"\n")
  postcov <- cov(mcmc.out$posterior)
  colnames(postcov) <- paste("theta",seq(1,mcmc.out$dim))
  rownames(postcov) <- paste("theta",seq(1,mcmc.out$dim))
  postsd <- as.table(postcov)
  print(postcov)
  cat("\n")
}

