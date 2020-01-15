RSLDiff<-function(data.one=NULL,data.two=NULL,modeldat.one=NULL,modeldat.two=NULL,num.pred=1000)
{
  dir.create(paste0("results/","RSLDifferences"),showWarnings = F)
  dataname.one<-data.one$dataname
  dataname.two<-data.two$dataname
  load(paste0("modeloutput/",dataname.one,"/run.mcmc.RData"))
  mcmcmodel.one <- run.mcmc$BUGSoutput$sims.list
  load(paste0("modeloutput/",dataname.two,"/run.mcmc.RData"))
  mcmcmodel.two <- run.mcmc$BUGSoutput$sims.list
  
  #Use the predictive distibution of the Gaussian process to determine the maximum Z score
  #This needs to be carried out on a grid of x values.
  N.grid.one <- modeldat.one$Ngrid
  x.grid.one <- modeldat.one$xstar
  xstar.one<-modeldat.one$xstar
  
  N.grid.two <- modeldat.two$Ngrid
  x.grid.two <- modeldat.two$xstar
  xstar.two<-modeldat.two$xstar
  
  xpred.one<-round(modeldat.one$year.grid,digits=0)
  xpred.two<-round(modeldat.two$year.grid,digits=0)
  
  index.year.one<-match(xpred.two,xpred.one)
  index.year.one<-index.year.one[which(!is.na(index.year.one))]
  
  index.year.two<-match(xpred.one,xpred.two)
  index.year.two<-index.year.two[which(!is.na(index.year.two))]
  
  
  #Set up the matrix that will contain the predictions
  pred.one.temp<-matrix(NA,ncol=N.grid.one,nrow=num.pred)
  pred.two.temp<-matrix(NA,ncol=N.grid.two,nrow=num.pred)
  
  K.gw.one<-matrix(NA, N.grid.one, N.grid.one)
  K.gw.two<-matrix(NA, N.grid.two, N.grid.two)
  
  ########Initialize quadrature for the integration########
  L=30    ## this sets the precision of the integration quadrature (higher is better but more computationally expensive)
  index=1:L        
  cosfunc=cos(((2*index-1)*pi)/(2*L))
  
  quad1.one=array(dim=c(nrow=N.grid.one,ncol=N.grid.one,L))
  quad2.one=array(dim=c(nrow=N.grid.one,ncol=N.grid.one,L))
  quad1.two=array(dim=c(nrow=N.grid.two,ncol=N.grid.two,L))
  quad2.two=array(dim=c(nrow=N.grid.two,ncol=N.grid.two,L))
  
  for(j in 1:N.grid.one)
  {   for(k in 1:N.grid.one) 
  { 
  quad1.one[k,j,]=abs((x.grid.one[k]*cosfunc/2)+(x.grid.one[k]/2)-xstar.one[j])^1.99
  quad2.one[k,j,]=((x.grid.one[k]/2)*(pi/L))*(sqrt(1-cosfunc^2))
  }
  }
  for(j in 1:N.grid.two)
  {   for(k in 1:N.grid.two) 
  {
  quad1.two[k,j,]=abs((x.grid.two[k]*cosfunc/2)+(x.grid.two[k]/2)-xstar.two[j])^1.99
  quad2.two[k,j,]=((x.grid.two[k]/2)*(pi/L))*(sqrt(1-cosfunc^2))
  
  }
  }
  
  for(i in 1:N.grid.one) { 
    for(j in 1:N.grid.one) {
      K.gw.one[j,i]<-sum((mean(mcmcmodel.one$p)^quad1.one[j,i,])*quad2.one[j,i,])  #### Quadrature function 
    }
  }
      for(i in 1:N.grid.two) { 
        for(j in 1:N.grid.two) {
          
      K.gw.two[j,i]<-sum((mean(mcmcmodel.two$p)^quad1.two[j,i,])*quad2.two[j,i,])  #### Quadrature function 
      
    } #End j loop 
  } #End i loop  
  
  pred.one<-pred.two<-pred.diff <- matrix(NA,ncol=length(index.year.one),nrow=num.pred)
  
  #Sampling from the predictive distibution
  for(i in 1:num.pred) {
    pred.one.temp[i,] <- mcmcmodel.one$beta0[i]+K.gw.one%*%mcmcmodel.one$K.w.inv[i,,]%*%mcmcmodel.one$w.m[i,]
    pred.two.temp[i,] <- mcmcmodel.two$beta0[i]+K.gw.two%*%mcmcmodel.two$K.w.inv[i,,]%*%mcmcmodel.two$w.m[i,]
    
    pred.one[i,]<-pred.one.temp[i,index.year.one]
    pred.two[i,]<-pred.two.temp[i,index.year.two]
    pred.diff[i,]<-pred.one[i,]-pred.two[i,]
  }
  
  diff.mean<-apply(pred.diff,2,mean)
  pred.one.mean<-apply(pred.one,2,mean)
  pred.two.mean<-apply(pred.two,2,mean)
  
  diff.u95<-apply(pred.diff,2,quantile,prob=0.975,na.rm=T)
  diff.l95<-apply(pred.diff,2,quantile,prob=0.025,na.rm=T)
  
  
  pdf(paste0("results/","RSLDifferences/","Differences",".pdf", sep = ""),height=6,width=10)
  
  plot(xpred.one[index.year.one],diff.mean,col="red",type="l",ylab="RSL Differences",xlab="Year CE ", main=" ",ylim=range(c(diff.l95,diff.u95)),xlim=range(xpred.one[index.year.one]),cex.axis=1.1)
  polygon(c(xpred.one[index.year.one], rev(xpred.one[index.year.one])),c(diff.l95,rev(diff.u95)), col=rgb(0.5,0,0.5,alpha=0.3),border="white")
  #polygon(rate_xxx, rate_zzz, col=colors()[596], border =colors()[596])
  abline(h=0,lty=2)
  dev.off()
  
   pdf(paste0("results/","RSLDifferences/","RSLDifferences",".pdf", sep = ""))
   
   plot(-10,-10,xlim=c(min(data.one$KYear*1000,data.two$KYear*1000),max(data.one$KYear*1000,data.two$KYear*1000)),ylim=c(min(data.one$RSL,data.two$RSL),max(data.one$RSL,data.two$RSL)),ylab="Relative Sea Level (m)",col=colors()[31],xlab="Year CE",main=paste0("RSL Differences (", modeldat.one$interval," ", "yr time steps)"),cex.axis=0.9)
   points(data.one$KYear*1000,data.one$RSL,col=rgb(1,0,0,alpha=0.5),lwd=0.5)
   points(data.two$KYear*1000,data.two$RSL,col=rgb(0,0,1,alpha=0.5),lwd=0.5)
   polygon(c(xpred.one[index.year.one],rev(xpred.one[index.year.one])),c(pred.one.mean,rev(pred.two.mean)),col=rgb(0.5,0,0.5,alpha=0.3),border="white")
   lines(xpred.one[index.year.one],pred.one.mean,col="red",lwd=1.5)
   lines(xpred.one[index.year.one],pred.two.mean,col="blue",lwd=1.5)
   legend("topleft",legend=c(dataname.one,dataname.two),col=c("red","blue"),lty=c(1,1),bty="n")
   dev.off()
  
  
  rsldiff<-cbind(year.pred=xpred.one[index.year.one],diff.mean,diff.l95,diff.u95)
  write.csv(rsldiff,file="results/RSLDifferences/Differences.csv")
  
  
  dir.create(paste0("modeloutput/","RSLDifferences"),showWarnings = F)
  pred.differences<-list(pred.one=pred.one.temp,pred.two=pred.two.temp,pred.diff=pred.diff)
  save(pred.differences, file=paste0("modeloutput/","RSLDifferences","/pred.differences.rda"))
  cat("Posteriors samples for RSL predictions and differences saved to modeloutput folder", "\n")
}

Diffresults<-function(data.one=NULL,data.two=NULL,modeldat.one=NULL,modeldat.two=NULL)
{
  load(paste0("modeloutput/","RSLDifferences","/pred.differences.rda"))

  dataname.one<-data.one$dataname
  dataname.two<-data.two$dataname
  yfit.one=pred.differences$pred.one
  xpred.one<-modeldat.one$year.grid
  meanfit.one<-apply(yfit.one,2,mean)
  u95.one<-meanfit.one+(2*apply(yfit.one,2,sd))
  l95.one<-meanfit.one-(2*apply(yfit.one,2,sd))
  u68.one<-meanfit.one+(1*apply(yfit.one,2,sd))
  l68.one<-meanfit.one-(1*apply(yfit.one,2,sd))
  
  yfit.two=pred.differences$pred.two
  xpred.two<-modeldat.two$year.grid
  meanfit.two<-apply(yfit.two,2,mean)
  u95.two<-meanfit.two+(2*apply(yfit.two,2,sd))
  l95.two<-meanfit.two-(2*apply(yfit.two,2,sd))
  u68.two<-meanfit.two+(1*apply(yfit.two,2,sd))
  l68.two<-meanfit.two-(1*apply(yfit.two,2,sd))
  
  dir.create(paste0("fig/","RSLDifferences"),showWarnings = F)
  
  pdf(paste0("fig/","RSLDifferences/","RSLDifferences",".pdf", sep = ""))
  
  plot(-10,-10,xlim=c(min(data.one$KYear*1000,data.two$KYear*1000),max(data.one$KYear*1000,data.two$KYear*1000)),ylim=c(min(data.one$RSL,data.two$RSL),max(data.one$RSL,data.two$RSL)),ylab="Relative Sea Level (m)",col=colors()[31],xlab="Year CE",main=paste0("RSL Differences(", modeldat.one$interval," ", "yr time steps)"),cex.axis=0.9)
  points(data.one$KYear*1000,data.one$RSL,col=rgb(1,0,0,alpha=0.5),lwd=0.5)
  points(data.two$KYear*1000,data.two$RSL,col=rgb(0,0,1,alpha=0.5),lwd=0.5)
  #polygon(c(xpred.one,rev(xpred.one)),c(meanfit.one,rev(meanfit.two)),col=rgb(0.5,0,0.5,alpha=0.3),border="white")
  lines(xpred.one,meanfit.one,col="red",lwd=1.5)
  lines(xpred.two,meanfit.two,col="blue",lwd=1.5)
  legend("topleft",legend=c(dataname.one,dataname.two),col=c("red","blue"),lty=c(1,1),bty="n")
  dev.off()
  
  cat(paste0("Plot of RSL differences between ",dataname.one," and ",dataname.two, " saved to fig folder"),sep="\n")
  
  dir.create(paste0("results/","RSLDifferences"),showWarnings = F)
  
  diff.mean<-apply(pred.differences$pred.diff,2,mean)
  diff.u95<-apply(pred.differences$pred.diff,2,quantile,prob=0.975,na.rm=T)
  diff.l95<-apply(pred.differences$pred.diff,2,quantile,prob=0.025,na.rm=T)
  
  modelpredone<-cbind(year.pred=xpred.one,RSL.pred=meanfit.one,l95=l95.one,u95=u95.one,l68=l68.one,u68=u68.one)
  write.csv(modelpredone,file=paste0("results/RSLDifferences/",dataname.one,"Predictions.csv"))
  
  modelpredtwo<-cbind(year.pred=xpred.two,RSL.pred=meanfit.two,l95=l95.two,u95=u95.two,l68=l68.two,u68=u68.two)
  write.csv(modelpredtwo,file=paste0("results/RSLDifferences/",dataname.two,"Predictions.csv"))
  
  rsldiff<-cbind(year.pred=xpred.two,diff.mean,diff.l95,diff.u95)
  write.csv(rsldiff,file="results/RSLDifferences/Differences.csv")
  
  cat(paste0("Spredsheets of RSL predictions and differences between ",dataname.one," and ",dataname.two, " saved to results folder"),sep="\n")
}


SLDiff<-function(data.one=NULL,data.two=NULL,modeldat.one=NULL,modeldat.two=NULL,num.pred=1000)
{
  dataname.one<-data.one$dataname
  dataname.two<-data.two$dataname
  load(paste0("modeloutput/",dataname.one,"/run.mcmc.RData"))
  mcmcmodel.one <- run.mcmc$BUGSoutput$sims.list
  load(paste0("modeloutput/",dataname.two,"/run.mcmc.RData"))
  mcmcmodel.two <- run.mcmc$BUGSoutput$sims.list
  
  #Use the predictive distibution of the Gaussian process to determine the maximum Z score
  #This needs to be carried out on a grid of x values.
  N.grid <- modeldat.one$Ngrid
  x.grid.one <- modeldat.one$xstar
  xstar.one<-modeldat.one$xstar
  
  x.grid.two <- modeldat.two$xstar
  xstar.two<-modeldat.two$xstar
  
  #Set up the matrix that will contain the predictions
  pred.one<-pred.two<-pred.diff <- matrix(NA,ncol=N.grid,nrow=num.pred)
  K.gw.one<-K.gw.two<-matrix(NA, N.grid, N.grid)
  
  ########Initialize quadrature for the integration########
  L=30    ## this sets the precision of the integration quadrature (higher is better but more computationally expensive)
  index=1:L        
  cosfunc=cos(((2*index-1)*pi)/(2*L))
  
  quad1.one=array(dim=c(nrow=N.grid,ncol=N.grid,L))
  quad2.one=array(dim=c(nrow=N.grid,ncol=N.grid,L))
  quad1.two=array(dim=c(nrow=N.grid,ncol=N.grid,L))
  quad2.two=array(dim=c(nrow=N.grid,ncol=N.grid,L))
  
  for(j in 1:N.grid)
  {   for(k in 1:N.grid) 
  { 
    quad1.one[k,j,]=abs((x.grid.one[k]*cosfunc/2)+(x.grid.one[k]/2)-xstar.one[j])^1.99
    quad2.one[k,j,]=((x.grid.one[k]/2)*(pi/L))*(sqrt(1-cosfunc^2))
    
    quad1.two[k,j,]=abs((x.grid.two[k]*cosfunc/2)+(x.grid.two[k]/2)-xstar.two[j])^1.99
    quad2.two[k,j,]=((x.grid.two[k]/2)*(pi/L))*(sqrt(1-cosfunc^2))
    
  }
  }
  
  for(i in 1:N.grid) { 
    for(j in 1:N.grid) {
      K.gw.one[j,i]<-sum((mean(mcmcmodel.one$p)^quad1.one[j,i,])*quad2.one[j,i,])  #### Quadrature function 
      K.gw.two[j,i]<-sum((mean(mcmcmodel.two$p)^quad1.two[j,i,])*quad2.two[j,i,])  #### Quadrature function 
      
    } #End j loop 
  } #End i loop  
  
  
  #Sampling from the predictive distibution
  for(i in 1:num.pred) {
    pred.one[i,] <- mcmcmodel.one$beta0[i]+K.gw.one%*%mcmcmodel.one$K.w.inv[i,,]%*%mcmcmodel.one$w.m[i,]
    pred.two[i,] <- mcmcmodel.two$beta0[i]+K.gw.two%*%mcmcmodel.two$K.w.inv[i,,]%*%mcmcmodel.two$w.m[i,]
    pred.diff[i,]<-pred.one[i,]-pred.two[i,]
  }
  
  dir.create(paste0("modeloutput/","SLDifferences"),showWarnings = F)
  pred.differences<-list(pred.one=pred.one,pred.two=pred.two,pred.diff=pred.diff)
  save(pred.differences, file=paste0("modeloutput/","SLDifferences","/pred.differences.rda"))
  cat("Posteriors samples for SL predictions and differences saved to modeloutput folder", "\n")
}

Diffresults.SL<-function(data.one=NULL,data.two=NULL,modeldat.one=NULL,modeldat.two=NULL)
{
  load(paste0("modeloutput/","SLDifferences","/pred.differences.rda"))
  
  dataname.one<-data.one$dataname
  dataname.two<-data.two$dataname
  yfit.one=pred.differences$pred.one
  xpred.one<-modeldat.one$year.grid
  meanfit.one<-apply(yfit.one,2,mean)
  u95.one<-meanfit.one+(2*apply(yfit.one,2,sd))
  l95.one<-meanfit.one-(2*apply(yfit.one,2,sd))
  u68.one<-meanfit.one+(1*apply(yfit.one,2,sd))
  l68.one<-meanfit.one-(1*apply(yfit.one,2,sd))
  
  yfit.two=pred.differences$pred.two
  xpred.two<-modeldat.two$year.grid
  meanfit.two<-apply(yfit.two,2,mean)
  u95.two<-meanfit.two+(2*apply(yfit.two,2,sd))
  l95.two<-meanfit.two-(2*apply(yfit.two,2,sd))
  u68.two<-meanfit.two+(1*apply(yfit.two,2,sd))
  l68.two<-meanfit.two-(1*apply(yfit.two,2,sd))
  
  dir.create(paste0("fig/","SLDifferences"),showWarnings = F)
  
  pdf(paste0("fig/","SLDifferences/","SLDifferences",".pdf", sep = ""))
  
  plot(-10,-10,xlim=c(min(data.one$KYear*1000,data.two$KYear*1000),max(data.one$KYear*1000,data.two$KYear*1000)),ylim=c(min(data.one$RSL,data.two$RSL),max(data.one$RSL,data.two$RSL)),ylab="Relative Sea Level (m)",col=colors()[31],xlab="Year CE",main=paste0("SL Differences(", modeldat.one$interval," ", "yr time steps)"),cex.axis=0.9)
  points(data.one$KYear*1000,data.one$RSL,col=rgb(1,0,0,alpha=0.5),lwd=0.5)
  points(data.two$KYear*1000,data.two$RSL,col=rgb(0,0,1,alpha=0.5),lwd=0.5)
  polygon(c(xpred.one,rev(xpred.one)),c(meanfit.one,rev(meanfit.two)),col=rgb(0.5,0,0.5,alpha=0.3),border="white")
  lines(xpred.one,meanfit.one,col="red",lwd=1.5)
  lines(xpred.two,meanfit.two,col="blue",lwd=1.5)
  legend("topleft",legend=c(dataname.one,dataname.two),col=c("red","blue"),lty=c(1,1),bty="n")
  dev.off()
  
  cat(paste0("Plot of SL differences between ",dataname.one," and ",dataname.two, " saved to fig folder"),sep="\n")
  
  dir.create(paste0("results/","SLDifferences"),showWarnings = F)
  
  diff.mean<-apply(pred.differences$pred.diff,2,mean)
  diff.u95<-apply(pred.differences$pred.diff,2,quantile,prob=0.975,na.rm=T)
  diff.l95<-apply(pred.differences$pred.diff,2,quantile,prob=0.025,na.rm=T)
  
  modelpredone<-cbind(year.pred=xpred.one,SL.pred=meanfit.one,l95=l95.one,u95=u95.one,l68=l68.one,u68=u68.one)
  write.csv(modelpredone,file=paste0("results/SLDifferences/",dataname.one,"Predictions.csv"))
  
  modelpredtwo<-cbind(year.pred=xpred.two,SL.pred=meanfit.two,l95=l95.two,u95=u95.two,l68=l68.two,u68=u68.two)
  write.csv(modelpredtwo,file=paste0("results/SLDifferences/",dataname.two,"Predictions.csv"))
  
  SLdiff<-cbind(year.pred=xpred.two,diff.mean,diff.l95,diff.u95)
  write.csv(SLdiff,file="results/SLDifferences/Differences.csv")
  
  cat(paste0("Spredsheets of SL predictions and differences between ",dataname.one," and ",dataname.two, " saved to results folder"),sep="\n")
}


RateDiff<-function(data.one=NULL,
                   data.two=NULL,
                   modeldat.one=NULL,
                   modeldat.two=NULL,
                   num.pred=1000,
                   main.diff="Rate Difference",
                   main.rate="Rates Overlayed")
{
  dataname.one<-data.one$dataname
  dataname.two<-data.two$dataname
  load(paste0("modeloutput/",dataname.one,"/run.mcmc.RData"))
  mcmcmodel.one <- run.mcmc$BUGSoutput$sims.list
  load(paste0("modeloutput/",dataname.two,"/run.mcmc.RData"))
  mcmcmodel.two <- run.mcmc$BUGSoutput$sims.list
  
  #Use the predictive distibution of the Gaussian process to determine the maximum Z score
  #This needs to be carried out on a grid of x values.
  N.grid <- modeldat.one$Ngrid
  x.grid.one <- modeldat.one$xstar
  xstar.one<-modeldat.one$xstar
  
  x.grid.two <- modeldat.two$xstar
  xstar.two<-modeldat.two$xstar
  
  xpred.one<-round(modeldat.one$year.grid,digits=0)
  xpred.two<-round(modeldat.two$year.grid,digits=0)
  
  index.year.one<-match(xpred.two,xpred.one)
  index.year.one<-index.year.one[which(!is.na(index.year.one))]
  
  index.year.two<-match(xpred.one,xpred.two)
  index.year.two<-index.year.two[which(!is.na(index.year.two))]
  
  #Set up the matrix that will contain the predictions
  pred.one<-pred.two<-pred.diff <- matrix(NA,ncol=length(index.year.one),nrow=num.pred)
  
  
  #Sampling from the predictive distibution
  for(i in 1:num.pred) {
    pred.one[i,] <- mcmcmodel.one$w.m[i,index.year.one]
    pred.two[i,] <- mcmcmodel.two$w.m[i,index.year.two]
    pred.diff[i,]<-pred.one[i,]-pred.two[i,]
  }
  
  dir.create(paste0("modeloutput/","RateDifferences"),showWarnings = F)
  rate.differences<-list(pred.one=pred.one,pred.two=pred.two,pred.diff=pred.diff)
  save(rate.differences, file=paste0("modeloutput/","RateDifferences","/rate.differences.rda"))
  cat("Posteriors samples for rates and differences saved to modeloutput folder", "\n")
  
  dir.create(paste0("results/","RateDifferences"),showWarnings = F)
  
  diff.mean<-apply(rate.differences$pred.diff,2,mean)
  diff.u95<-apply(rate.differences$pred.diff,2,quantile,prob=0.975,na.rm=T)
  diff.l95<-apply(rate.differences$pred.diff,2,quantile,prob=0.025,na.rm=T)
  
  Ratediff<-cbind(year.pred=xpred.one[index.year.one],diff.mean,diff.l95,diff.u95)
  write.csv(Ratediff,file="results/RateDifferences/Differences.csv")
  
  cat(paste0("Spredsheets of SL predictions and differences between ",dataname.one," and ",dataname.two, " saved to results folder"),sep="\n")
  
  rates.one=rate.differences$pred.one
  mean.one<-apply(rates.one,2,mean)
  u95.one<-apply(rates.one,2,quantile,prob=0.975,na.rm=T)
  l95.one<-apply(rates.one,2,quantile,prob=0.025,na.rm=T)

  rates.two=rate.differences$pred.two
  mean.two<-apply(rates.two,2,mean)
  u95.two<-apply(rates.two,2,quantile,prob=0.975,na.rm=T)
  l95.two<-apply(rates.two,2,quantile,prob=0.025,na.rm=T)
  

  pdf(paste0("results/","RateDifferences/","RateDifferences",".pdf", sep = ""),height=6,width=10)
  
  plot(xpred.one[index.year.one],mean.one,col="red",type="l",ylab="Rate (mm/yr)",xlab="Year AD ", main=main.rate,ylim=c(min(l95.one, l95.two),max(u95.one,u95.two)),xlim=range(xpred.one[index.year.one]),cex.axis=1.1)
  lines(xpred.two[index.year.two],mean.two,col="blue")
  #polygon(c(xpred.one,rev(xpred.one)),c(mean.one,rev(mean.two)),col=rgb(0.5,0,0.5,alpha=0.3),border="white")
  
  polygon(c(xpred.one[index.year.one], rev(xpred.one[index.year.one])),c(l95.one,rev(u95.one)), col=rgb(1,0,0,alpha=0.3),border="white")
  polygon(c(xpred.two[index.year.two], rev(xpred.two[index.year.two])),c(l95.two,rev(u95.two)), col=rgb(0,0,1,alpha=0.3),border="white")
  #abline(h=0)
  #polygon(rate_xxx, rate_zzz, col=colors()[596], border =colors()[596])
  dev.off()
  
  pdf(paste0("results/","RateDifferences/","Differences",".pdf", sep = ""),height=6,width=10)
  
  plot(xpred.one[index.year.one],diff.mean,col="red",type="l",ylab="Rate Differences",xlab="Year AD ", main=main.diff,ylim=range(c(diff.l95,diff.u95)),xlim=range(xpred.one[index.year.one]),cex.axis=1.1)
  polygon(c(xpred.one[index.year.one], rev(xpred.one[index.year.one])),c(diff.l95,rev(diff.u95)), col=rgb(0.5,0,0.5,alpha=0.3),border="white")
  #polygon(rate_xxx, rate_zzz, col=colors()[596], border =colors()[596])
  abline(h=0,lty=2)
  dev.off()
  
 }
