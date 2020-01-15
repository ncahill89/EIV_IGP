#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
ConstructMCMCArray <- function(# Read in JAGS objects
  ###  Read in JAGS objects and constructs \code{mcmc.array}, 
  ### which is saved to \code{output.dir}.
  ### This function can only be run after finising the mcmc sampling (after function \code{\link{RunMCMC}} has completed).
  ChainIDs = ChainNums, ##<< Optional: specify which chains to include
  ## (to use when you want to exclude a chain that crashed, or which has not finished yet).
  n.samplestot.max = 15000, ##<< Maximum number of posterior samples to save
  output.dir = NULL, ##<< Directory where MCMC output was stored and will be stored. 
  core.run=FALSE,
  data.raw
){
  
  dataname<-data.raw$dataname
  
  if (is.null(output.dir)){
    output.dir<-file.path(paste0("modeloutput/",dataname))
  }
  
  # now combine the JAGS files into one mcmc.array
  n.chains <- length(ChainIDs)
  if (n.chains==1){
    cat("You need at least two chains!\n")
    return()
  }
  jags.dir <- file.path(output.dir, "temp.JAGSobjects/")
  
  cat("Reading in JAGS output files from", jags.dir, "\n")
  chain  <- ifelse(length(ChainIDs)==1,ChainIDs,ChainIDs[1])
  load(file.path(jags.dir, paste0("jags_mod", chain, ".Rdata"))) 
  n.sim <- dim(mod.upd$BUGSoutput$sims.array)[1]
  n.par <- dim(mod.upd$BUGSoutput$sims.array)[3]
  mcmc.array <- array(NA, c(n.sim, n.chains, n.par))
  dimnames(mcmc.array) <- list(NULL, NULL, names(mod.upd$BUGSoutput$sims.array[1,1,]))
  for (chain in 1:n.chains){
    chain_saved <- ifelse(length(ChainIDs)==1,1,ChainIDs[chain])
    cat(paste("Reading in chain number ", chain_saved, sep = ""), "\n")
    load(file.path(output.dir, "temp.JAGSobjects", paste0("jags_mod", chain_saved, ".Rdata"))) 
    mcmc.array[1:n.sim,chain, ] = mod.upd$BUGSoutput$sims.array[,1,]
    
  }

  if (n.sim > n.samplestot.max){
    mcmc.array <- mcmc.array[seq(1, n.sample.max, length.out = n.sample.max), , ]  
  }
  
  if(!core.run)
  save(mcmc.array, file = file.path(output.dir, paste0("mcmc.array", ".rda"))) 
  if(core.run)
  save(mcmc.array, file = file.path(output.dir, paste0("mcmc.array.core", ".rda"))) 
  
  cat("mcmc.array saved to", output.dir, "\n")
  
  

  return(invisible())
}
#----------------------------------------------------------------------   
# The End!
