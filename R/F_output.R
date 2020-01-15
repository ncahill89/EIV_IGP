IGPests<-function(data.raw=NULL,
                  interval = 25)
{
  dataname<-data.raw$dataname
  
  # Get model data
  modeldat <- IGPdata(data.raw = data.raw,
                      interval = interval)
  # Get model output
  load(paste0("modeloutput/",dataname,"/mcmc.array.rda"))
  n_iter <- length(mcmc.array[,1,paste0("beta0")])
  
  #Get predictions on a grid of x values.
  N.grid <- modeldat$Ngrid
  x.grid <- modeldat$xstar
  xstar <- modeldat$xstar
  Dist <- modeldat$Dist
  
  #Set up the matrix that will contain the estimates
  pred <- matrix(NA,ncol=N.grid,nrow=n_iter)
  K.gw<-K<-K.w.inv<-array(NA,c(n_iter, N.grid, N.grid))
  
  ########Initialize quadrature for the integration########
  L=30    ## this sets the precision of the integration quadrature (higher is better but more computationally expensive)
  index=1:L        
  cosfunc=cos(((2*index-1)*pi)/(2*L))
  
  quad1=array(dim=c(nrow=N.grid,ncol=N.grid,L))
  quad2=array(dim=c(nrow=N.grid,ncol=N.grid,L))
  
  for(j in 1:N.grid)
  {   
    for(k in 1:N.grid) 
  { 
      quad1[k,j,]=abs((x.grid[k]*cosfunc/2)+(x.grid[k]/2)-xstar[j])^1.99
      quad2[k,j,]=((x.grid[k]/2)*(pi/L))*(sqrt(1-cosfunc^2))
  }
  }
  
  
  #Get posterior samples of rates
  w.ms<-array(NA, c(n_iter,N.grid))
  for(j in 1:N.grid)
  {
    w.ms[,j]<-mcmc.array[,1,paste0("w.m[",j,"]")]
  }
  
  #Get estimates
  for(i in 1:n_iter) {
    for(k in 1:N.grid) { 
      for(j in 1:N.grid) {
        K.gw[i,j,k]<-sum((mcmc.array[i,1,paste0("p")]^quad1[j,k,])*quad2[j,k,])  #### Quadrature function 
      } #End j loop 
    } #End k loop  
    
    K[i,,]<-mcmc.array[i,1,paste0("p")]^(Dist^1.99)
    K.w.inv[i,,]<-solve(K[i,,])
    pred[i,] <- mcmc.array[i,1,paste0("beta0")]+K.gw[i,,]%*%K.w.inv[i,,]%*%w.ms[i,]
  } #End i loop
  
  dydt<-array(NA,c(n_iter,N.grid))
  mean.dydt<-rep(NA,n_iter)
  for(i in 1:n_iter){
    dydt[i,] <- w.ms[i,]
    mean.dydt[i]<-mean(dydt[i,])
  }
  
  mean.rate <- mean(mean.dydt)
  sd.rate <- sd(mean.dydt)
  u95.rate <- quantile(mean.dydt,probs=0.975)
  l95.rate <- quantile(mean.dydt,probs=0.025)
  
  EstsandRates<-list(pred=pred,dydt=dydt,mean.rate=mean.rate,sd.rate=sd.rate,l95.rate=l95.rate,u95.rate=u95.rate)
  save(EstsandRates, file=paste0("modeloutput/",dataname,"/EstsandRates.rda"))
  cat("EIV-IGP posterior samples for estimates and rates saved to modeloutput folder", "\n")
}