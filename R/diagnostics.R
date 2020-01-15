gr_diag<-function(mcmc.array,
                  pars.check=c("p","sigma.g","beta0")){
R <- rep(NA,length(pars.check))
p<-0
for (parname in pars.check){
  p <- p+1 # index
  mcmc.array.temp <- mcmc.array[,,parname]
  mcmc <- mcmc.list()
  
for (chain in 1:dim(mcmc.array.temp)[2]){
  mcmc[[chain]] <- as.mcmc(mcmc.array.temp[,chain])
}
  r<- gelman.diag(mcmc, autoburnin = FALSE, transform = F)$psrf
  R[p] <-r[,"Point est."]
  
}


names(R) <- pars.check


if (length(R[R>1.1])>0){
  cat(paste("Poor/no convergence for:", names(R[R>1.1]), "(R = ", round(R[R>1.1],3), ")", "\n"))
}
else
  cat(paste0("Rhat looks good, no  convergence issues indicated for checked parameters \n"))

if(length(R[R>1.1])>0)
  return(-1)
else
  return(0)

}


eff_size<-function(mcmc.array,
                   pars.check=c("p","sigma.g","beta0")){
  ESS <- rep(NA,length(pars.check))
  p<-0
  for (parname in pars.check){
    p <- p+1 # index
    mcmc.array.temp <- mcmc.array[,,parname]
    mcmc <- mcmc.list()
    
    for (chain in 1:dim(mcmc.array.temp)[2]){
      mcmc[[chain]] <- as.mcmc(mcmc.array.temp[,chain])
    }
    es<- effectiveSize(mcmc)
    ESS[p] <-es/(dim(mcmc.array)[1]*dim(mcmc.array)[2])
    
  }
  names(ESS) <- pars.check
  if (length(ESS[ESS<0.1])>0){
    cat(paste0("Effective sample size is less than 10% of total iterations for parameter:"," ", names(ESS[ESS<0.1])," ", "(",round(ESS[ESS<0.1],3)*100,"%",")", "\n"))
    cat(paste0("Additional thinning may be required! \n"))
    }
  else
    cat(paste0("No apparent autocorrelation issues for checked parameters. \n"))

  
}

mcse<-function(mcmc.array,
               pars.check=c("p","sigma.g","beta0")){
  MCSE <- rep(NA,length(pars.check))
  p<-0
  for (parname in pars.check){
    p <- p+1 # index
    mcmc.array.temp <- mcmc.array[,,parname]
    mcmc <- mcmc.list()
    
    for (chain in 1:dim(mcmc.array.temp)[2]){
      mcmc[[chain]] <- as.mcmc(mcmc.array.temp[,chain])
    }
    es<- effectiveSize(mcmc)
    MCSE[p] <-(sd(mcmc.array.temp)/es)/sd(mcmc.array.temp)
    
  }
  names(MCSE) <- pars.check
  if (length(MCSE[MCSE>0.1])>0){
    cat(paste0("The Monte Carlo standard error is greater than 10% of the posterior standard deviation for parameter:"," ", names(MCSE[MCSE>0.1])," ", "(",round(MCSE[MCSE>0.1],3)*100,"%",")", "\n"))
    cat(paste0("Sampling error variation appears too large! \n"))
  }
  else
    cat(paste0("The accuracy of the parameter estimation is adequate. \n"))
}




