RunIGPModel<-function(data.raw=NULL,
                      cor.p=0.2,
                      n.iter=30000,
                      n.burnin=10000,
                      n.thin=10,
                      ChainNums=seq(1,2),
                      fast=FALSE,
                      run.on.server=FALSE,
                      interval = 25){
  # Create a directory "modeloutput" in current working directory
  dataname<-data.raw$dataname
  dir.create("modeloutput", showWarnings = FALSE)
  dir.create(paste0("modeloutput/",dataname),showWarnings = F)
  output.dir<-file.path(paste0("modeloutput/",dataname))
  
  # Get model data
  modeldat <- IGPdata(data.raw = data.raw,
                      interval = interval)
  
  ###The necessary data
  jags.data <- list(n = modeldat$N,
                 m = modeldat$Ngrid,
                 P = modeldat$P,
                 D = modeldat$D,
                # L = modeldat$L,
                # ppi = modeldat$ppi,
                # cosfunc = modeldat$cosfunc,
                 Dist = modeldat$Dist,
                 xstar = modeldat$xstar,
                 quad1 = modeldat$quad1,
                 quad2 = modeldat$quad2,
                 kappa = 1.99,
                 cor.p = cor.p)    

  ###Paramaters to save 
  jags.pars <- c("beta0",
                 "sigma.g",
                 "p",
                 "w.m",
                 "mu.y",
                 "sigma.y",
                 "K.w.inv")

  ########Run the model########
  if(fast)
  {
    model.file="model/EIVIGPfast.txt"
  }
  
  if(!fast)
  {
    model.file="model/EIVIGP.txt"
  }
  
  
  if (run.on.server) {
    foreach(chainNum=ChainNums) %dopar% {
      cat(paste("Start chain ID ", chainNum), "\n")
      
      InternalRunOneChain(chainNum = chainNum,
                          jags.data = jags.data,
                          jags.pars = jags.pars,
                          n.burnin = n.burnin,
                          n.iter = n.iter,
                          n.thin = n.thin,
                          model.file = model.file,
                          output.dir = output.dir)
    } # end chainNums
  } else {
    for (chainNum in ChainNums){
      cat(paste("Start chain ID ", chainNum), "\n")
      
      InternalRunOneChain(chainNum = chainNum,
                          jags.data = jags.data,
                          jags.pars = jags.pars,
                          n.burnin = n.burnin,
                          n.iter = n.iter,
                          n.thin = n.thin,
                          model.file = model.file,
                          output.dir = output.dir)
      
    }
  }
  
  # contruct MCMC array
  ConstructMCMCArray(ChainIDs = ChainNums,data.raw = data.raw)
  
  # Get estimates
  IGPests(data.raw=data.raw,
          interval = interval)
  
} 


#-----------------------------------------------------
InternalRunOneChain <- function(#Do MCMC sampling
  ###Do MCMC sampling for one chain
  chainNum, ##<< Chain ID
  jags.data,
  jags.pars,
  n.burnin,
  n.iter,
  n.thin,
  output.dir,
  model.file
){
  # set seed before sampling the initial values
  set.seed.chain <- chainNum*209846
  # mcmc.info <- list(set.seed.chain = set.seed.chain, chainNum = chainNum)
  # mcmc.info.file <- file.path(mcmc.meta$general$output.dir, paste0("mcmc.info", filename.append, ".", chainNum, ".rda"))
  # save(mcmc.info, file = mcmc.info.file)
  dir.create(paste0(output.dir, "/temp.JAGSobjects/"),showWarnings=FALSE)
  jags.dir <- file.path(output.dir, "temp.JAGSobjects/")
  set.seed(set.seed.chain)
  temp <- rnorm(1)
  
  mod<-suppressWarnings(jags(data=jags.data,
            parameters.to.save=jags.pars,
            model.file=model.file,
            n.chains=1,
            n.iter=n.iter,
            n.burnin=n.burnin,
            n.thin=n.thin,
            DIC=FALSE,
            jags.seed = set.seed.chain))
  
  mod.upd <- mod
  save(mod.upd, file=file.path(output.dir, "temp.JAGSobjects", paste0("jags_mod", chainNum, ".Rdata")))
  cat(paste("MCMC results", " for chain ", chainNum, " written to folder temp.JAGSobjects in ", output.dir), "\n")
  
  
  cat(paste("Hooraah, Chain", chainNum, "has finished!"), "\n")
  return(invisible())
}
