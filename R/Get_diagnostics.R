

get_diagnostics<-function(data.raw,output.dir="modeloutput/")
{
  dataname<-data.raw$dataname
  load(paste0(output.dir,dataname,"/mcmc.array.rda"))
     
  pars.check=c("p","sigma.g","beta0")
  # Get gelman diagnostics (Rhat threshold = 1.1)
  # If gelman diagnostic fails then stop!
    gd<-gr_diag(mcmc.array,pars.check = pars.check)
    if(gd==-1)
    {
      cat("WARNING! Convergence issues, check trace plots \n")
    
    return(pars.check=pars.check)
    }
  # If gelman diagnostic passes then get other diagnostics  
    if(gd==0)
    {
    eff_size(mcmc.array,pars.check = pars.check)
    mcse(mcmc.array,pars.check = pars.check)
    
    return(list(pars.checked=pars.check))
    }
  
  
}
