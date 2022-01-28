source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
library(fields)
warm_cutoff = F
scale = "cen"
scaling <- list()

#be careful running this code, it might require up to 30 GB free memory
for(n in c(signal_tbb %>% filter(type=="model"))$signal){
  summary <- readRDS(paste0("processing/raw_data/", n, ".Rds"))

  if(n %in% signal_tbb$signal[which(grepl("highres",signal_tbb$signal))]){
    next
  }
  
  if(warm_cutoff){
    summary$temp <- summary$temp[,, index(summary$time[summary$time < warm_cutoff])]
    summary$time <- summary$time[summary$time < warm_cutoff]
  }
  print(n)
  
  dimnames(summary$temp) <- list(longitude=summary$lon,latitude=summary$lat,time=summary$time)
  #if(n!="IPSL"){#for IPSL it's too costly and memory will be crashed
  summary$temp <- summary$temp[c(length(summary$lon),1:length(summary$lon), 1),,]
  summary$lon <- c(max(summary$lon)-360,summary$lon, min(summary$lon)+360)
  #}
  min.res = get_min.res(tscale[[scale]])
  min.range = get_min.range(tscale[[scale]])
  max.hiat = get_max.hiat(tscale[[scale]])
  pages.meta <- readRDS("data/pages_meta_scaling.Rds")
  if((length(which(summary$lon < 0))>1)==FALSE){pages.meta <- pages.meta %>% mutate(Lon = case_when(Lon < 0 ~ Lon + 360, TRUE ~ Lon))}
   sim_data <- zoo(t(apply(summary$temp, 3, function(data_tmp){
      fields::interp.surface(list(x=summary$lon, y=summary$lat, z=data_tmp), data.frame(x=pages.meta$Lon, y=pages.meta$Lat))
      })), order.by=summary$time)

  rm(summary)
  gc()  
  
  res <- signal_tbb[which(signal_tbb$signal==n), ][["interp.res"]]
  time <- index(sim_data[,1])

  tmp <- apply(sim_data, 2, function(x) SpecMTM(MakeEquidistant(time, (coredata(ts(x))-mean(coredata(ts(x)))), dt=res), k=k, nw=nw, detrend=TRUE))
  fit <- lapply(tmp, function(x) SlopeFit(x, 1/max(tscale[[scale]]), 1/min(tscale[[scale]]), bDebug=T))
  scaling[[n]] <- lapply(fit, function(x) list(slope=x$slope, slopesd=x$slopesd))
  names(scaling[[n]]) <- pages.meta$Name  
}

tbb <- tibble()
for(n in names(scaling)){
  tbb <- rbind(tbb, list_to_tibble(scaling[[n]]) %>% unnest(data) %>% add_column(signal=n))
}
if(save){
  saveRDS(tbb, "processing/bilint_scalingcoeff.Rds")
}
