GetIGPRes<-function(data.raw=NULL,
                    interval = 25)
{

  dir.create("results", showWarnings = FALSE)
  dir.create(paste0("results/",data.raw$dataname),showWarnings = F)
  
  GIA <- data.raw$GIA
  BP <- data.raw$BP_age_scale
  
  load(paste0("modeloutput/",data.raw$dataname,"/EstsandRates.rda"))
  
  modeldat <- IGPdata(data.raw = data.raw,
                      interval = interval,
                      incl.errbounds = EstsandRates$incl.errbounds)
  

  
  pred_s <- suppressWarnings(as_tibble(EstsandRates$pred) %>% 
              rename_at(vars(everything()),~ as.character(modeldat$year.grid)) %>% 
              pivot_longer(everything(),names_to = "year") %>% 
              mutate(year = as.numeric(year)))
  SLestimates <- pred_s %>% 
                  group_by(year) %>% 
                  summarise(SL_est = mean(value),
                            SL_lwr = mean(value) - 2*(sd(value)),
                            SL_upr = mean(value) + 2*(sd(value))) %>% 
                  mutate(year = ifelse(rep(data.raw$BP_age_scale,length(SL_est)) == FALSE, year, 1950 - year))
  
  dydt_s <- suppressWarnings(as_tibble(EstsandRates$dydt) %>% 
            rename_at(vars(everything()),~ as.character(modeldat$year.grid)) %>% 
            pivot_longer(everything(),names_to = "year") %>% 
            mutate(year = as.numeric(year)))

  SLrates <- dydt_s %>% 
             group_by(year) %>% 
             summarise(rate_est = mean(value),
                       rate_lwr = mean(value) - 2*(sd(value)),
                       rate_upr = mean(value) + 2*(sd(value))) %>% 
             mutate(year = ifelse(rep(data.raw$BP_age_scale,length(rate_est)) == FALSE, year, 1950 - year))
  
  rate_mean <- mean(dydt_s$value)
  rate_mean_lwr <- quantile(dydt_s$value,probs = 0.025)
  rate_mean_upr <- quantile(dydt_s$value,probs = 0.975)
  
  write.csv(SLestimates,file=paste0("results/",data.raw$dataname,"/",ifelse(data.raw$GIA == FALSE,"RSL_Estimates.csv", "SL_Estimates.csv")))
  write.csv(SLrates,file=paste0("results/",data.raw$dataname,"/",ifelse(data.raw$GIA == FALSE, "RSL_Rates.csv","SL_Rates.csv")))
  cat(paste0("Spreadsheets containing ", ifelse(data.raw$GIA == FALSE, "RSL estimates ","GIA Corrected SL estimates "),"and rates for"," ",data.raw$dataname," ", "are saved in results folder","\n"))

  return(list(SLestimates = SLestimates,
              SLrates = SLrates, 
              modeldat = modeldat,
              rate_mean = rate_mean,
              rate_mean_lwr = rate_mean_lwr,
              rate_mean_upr = rate_mean_upr))
}

IGPResults<-function(data.raw=NULL,
                   interval = 25,
                   xlimits=NULL,
                   ylimits=NULL,
                   ratelimits=NULL)
{


  model_res <- GetIGPRes(data.raw=data.raw,
                         interval = interval)
  modeldat <- model_res$modeldat
  
  sl_dat <- model_res$SLestimates
  rate_dat <- model_res$SLrates
  data_to_plot <- data.raw$data_to_plot
  
  p1 <- ggplot(sl_dat, aes(x = year, y = SL_est))+
    geom_line() +
    geom_ribbon(aes(x = year, ymin = SL_lwr, ymax = SL_upr), alpha = 0.5) +
    geom_polygon(data = data_to_plot, aes(x = x, y = y,group = obs_index),alpha = 0.1) + 
    ylab(ifelse(data.raw$GIA == FALSE, "Relative Sea Level (m)", "Sea Level (m)")) +
    xlab(ifelse(data.raw$BP_age_scale == FALSE,"Year CE","Year BP")) + 
    ggtitle(ifelse(data.raw$GIA == FALSE,"RSL Estimates","SL Estimates")) + 
    theme_classic()
  if(data.raw$BP_age_scale == TRUE) {p1 <- p1 + scale_x_reverse()}
  
  suppressMessages(ggsave(paste0("fig/",data.raw$dataname,"/","Results_SL Estimates", ifelse(data.raw$GIA == FALSE,"","(GIA corrected)"),".pdf", sep = ""),p1, width = 7, height = 4))
  
  
  p2 <- ggplot(rate_dat, aes(x = year, y = rate_est))+
    geom_line() +
    geom_ribbon(aes(x = year, ymin = rate_lwr, ymax = rate_upr), alpha = 0.5) +
    ylab(ifelse(data.raw$GIA == FALSE, "Rate of RSL Change (mm/yr)", "Rate of SL Change (mm/yr)")) +
    xlab(ifelse(data.raw$BP_age_scale == FALSE,"Year CE","Year BP")) + 
    ggtitle(ifelse(data.raw$GIA == FALSE,"RSL Rates","SL Rates")) + 
    theme_classic()
  
  if(data.raw$BP_age_scale == TRUE) {p2 <- p2 + scale_x_reverse()}
  
  suppressMessages(ggsave(paste0("fig/",data.raw$dataname,"/","Results_Rate Estimates", ifelse(data.raw$GIA == FALSE,"","(GIA corrected)"),".pdf", sep = ""),p2, width = 7, height = 4))
  

  cat("Plots of estimates and rates saved in fig folder \n")
  
  return(list(mean_rate = model_res$rate_mean, 
              mean_rate_lwr = model_res$rate_mean_lwr,
              mean_rate_upr = model_res$rate_mean_upr))
}
