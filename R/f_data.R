dataprep<-function(data_path = NULL,
                   dataname="RSL Record",
                   GIA=FALSE,
                   rate.gia=NULL,
                   yocc=2010, 
                   BP_age_scale = FALSE)
{
  
  ####Create Directories
  dir.create("fig", showWarnings = FALSE)
  dir.create(paste0("fig/",dataname),showWarnings = F)

  ## read in data and make any necessary changes
  data <- read_csv(data_path,
                   col_types = cols())
  names(data)[grepl("AgeError",names(data))]<-"AgeError"
  names(data)[grepl("RSLError",names(data))]<-"RSLError"
  
  ########Set up the data########
  GIA <- rep(GIA, nrow(data))
  BP <- rep(BP_age_scale, nrow(data))
  
    data <- data %>% 
            mutate(x = ifelse(BP == FALSE, Age, 1950 - Age),
                   x_thousand = x/1000,
                   y = ifelse(GIA == FALSE,RSL,(((yocc/1000)-x_thousand)*rate.gia)+RSL),
                   var_x = (AgeError/1000)^2,
                   var_y = RSLError^2,
                   x_lwr = x - AgeError, 
                   x_upr = x + AgeError,
                   y_lwr = RSL - RSLError,
                   y_upr = RSL + RSLError,
                   y_1_lwr = ifelse(GIA == FALSE, y - RSLError,(((yocc/1000)-x_upr/1000)*rate.gia)+y_lwr),
                   y_2_lwr = ifelse(GIA == FALSE, y - RSLError,(((yocc/1000)-x_lwr/1000)*rate.gia)+y_lwr),
                   y_3_upr = ifelse(GIA == FALSE, y + RSLError,(((yocc/1000)-x_lwr/1000)*rate.gia)+y_upr),
                   y_4_upr = ifelse(GIA == FALSE, y + RSLError,(((yocc/1000)-x_upr/1000)*rate.gia)+y_upr),
                   x_1_upr = x + AgeError, 
                   x_2_lwr = x - AgeError,
                   x_3_lwr = x - AgeError, 
                   x_4_upr = x + AgeError)
            
  ########Setting up the covariance and precision matrices#######
  N <- nrow(data)
  V22 <- ifelse(GIA == FALSE,data$var_y,(((rate.gia^2)*data$var_x)+data$var_y)+data$var_x)
  V12<- ifelse(GIA == FALSE, 0, -rate.gia*data$var_x)
  V21<- ifelse(GIA == FALSE, 0, -rate.gia*data$var_x)
  V11 <- data %>% pull(var_x)
  V<-array(NA,c(2,2,nrow(data)))
  P<-array(NA,c(2,2,nrow(data)))
  for(i in 1:N)
  {
    V[,,i]<- matrix(c(V11[i],V12[i],V21[i],V22[i]),2,2)
    P[,,i]<-solve(V[,,i])
  }
  
  get_bounds <- data %>% 
            select(y_1_lwr:x_4_upr) %>% 
            mutate(obs_index = 1:n()) %>% 
            pivot_longer(cols = y_1_lwr:x_4_upr,
                         names_to = "bounds",
                         values_to = "value") %>% 
            mutate(bounds = replace(bounds, bounds %in% c("y_1_lwr","y_2_lwr","y_3_upr","y_4_upr"), "y"),
                   bounds = replace(bounds, bounds %in% c("x_1_upr","x_2_lwr","x_3_lwr","x_4_upr"), "x")) 

   x_bounds <- get_bounds %>% 
              filter(bounds == "x") 
 
   y_bounds <- get_bounds %>% 
                  filter(bounds == "y") 
 
  data_to_plot <- tibble(obs_index = x_bounds$obs_index, 
                           x = ifelse(rep(BP_age_scale,nrow(x_bounds)) == FALSE,x_bounds$value,1950 - x_bounds$value),
                           y = y_bounds$value)
  
  
  p <- ggplot(data_to_plot, aes(x = x, y = y))+
    geom_polygon(aes(group = obs_index),alpha = 0.3) + 
    geom_point(data = data, aes(x = Age, y = y), alpha = 0.6, pch = 1) +
    ylab(ifelse(GIA == FALSE, "Relative Sea Level (m)", "Sea Level (m)")) +
    xlab(ifelse(BP_age_scale == FALSE,"Year CE","Year BP")) + 
    ggtitle(ifelse(GIA == FALSE,"RSL Reconstruction","SL Reconstruction")) + 
    theme_classic()

    suppressMessages(ggsave(paste0("fig/",dataname,"/","Raw Data",ifelse(GIA[1] == FALSE,"","(GIA corrected)"),".pdf", sep = ""),p, width = 7, height = 4))
    cat("Plots of data saved to fig folder", "\n")

  
  return(list(data = data,
              data_to_plot = data_to_plot,
              P = P,
              dataname=dataname,
              GIA=GIA[1],
              BP_age_scale = BP_age_scale))

}


IGPdata<-function(data.raw = NULL,
                  year1 = NULL,
                  year2 = NULL,
                  interval = 25,
                  html.file = NULL,
                  incl.errbounds = TRUE,
                  upper = NULL,
                  lower = NULL)
{
  
  data <- data.raw$data
  
  ############# Set up the grid for the GP ###################
  nw=40      # Sets the min number of grid points for the derivative
  if(incl.errbounds){
  up <- max(data$x_upr)/1000
  low <- min(data$x_lwr)/1000
  xgrid=c(low,seq(min(data$x_thousand),max(data$x_thousand),by=(interval/1000)),up)
  }
  else{
    up = ifelse(is.null(upper),max(data$x_thousand),upper)
    low = ifelse(is.null(lower),min(data$x_thousand),lower)
    xgrid=c(low,seq(min(data$x_thousand),max(data$x_thousand),by=(interval/1000)),up)
  }
  Ngrid = length(xgrid)

  if(Ngrid<nw)
    stop("Grid length must be at least 40")
  
  else
    cat(paste0("Using a grid size of"," ",Ngrid," ", "and an interval width of"," ",interval," ","years \n"))
  
   ###Change data to lower zero for integration
   minx = min(data$x_thousand)
   x = data$x_thousand-minx
   xstar = xgrid - minx
   Dist <- rdist(xstar) ###Distance matrix required for the model 
   D <- cbind(x,data$y) ###Combine the x,y data for the model 
   
   ########Initialize quadrature for the integration########
   N <- nrow(data)
   L = 30    ## this sets the precision of the integration quadrature (higher is better but more computationally expensive)
   index=1:L        
   cosfunc=cos(((2*index-1)*pi)/(2*L))
   
   quad1=array(dim=c(nrow=N,ncol=Ngrid,L))
   quad2=array(dim=c(nrow=N,ncol=Ngrid,L))
   
   for(j in 1:Ngrid)
   {   for(k in 1:N) 
   { 
     quad1[k,j,]=abs((x[k]*cosfunc/2)+(x[k]/2)-xstar[j])^1.99
     quad2[k,j,]=((x[k]/2)*(pi/L))*(sqrt(1-cosfunc^2))
   }
   }

   P <- data.raw$P
   
   return(list(year.grid = xgrid*1000,
               xstar = xstar,
               N = N,
               Ngrid = Ngrid,
               D = D,
               P = P,
               Dist = Dist,
               quad1 = quad1,
               quad2 = quad2,
               cosfunc = cosfunc,
               ppi = pi,
               L = L,
               incl.errbounds=incl.errbounds,
               interval=interval,
               BP_age_scale = data.raw$BP_age_scale))
 }