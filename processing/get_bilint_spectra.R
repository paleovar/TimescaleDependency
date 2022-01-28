source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
warm_cutoff = F
scale = "cen"

for(n in c(signal_tbb %>% filter(type=="model"))$signal){
  summary <- readRDS(paste0("processing/raw_data/", n, ".Rds"))
  if(warm_cutoff){
  if(!n %in% c(signal_tbb %>% filter(type=="obs") %>% select(signal))$signal&!n %in% signal_tbb$signal[which(grepl("highres",signal_tbb$signal))])
  {
    summary$temp <- summary$temp[,, index(summary$time[summary$time < warm_cutoff])]
    summary$time <- summary$time[summary$time < warm_cutoff]
  }
  }
  
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
      interp.surface(list(x=summary$lon, y=summary$lat, z=data_tmp), data.frame(x=pages.meta$Lon, y=pages.meta$Lat))
      })), order.by=summary$time)
  #saveRDS(sim_data, paste0(pages.dir, "sim_data_", n, "_", scale, "_", round(min.res,3), "_", round(min.range,2), "_", round(max.hiat,2), ".Rds"))
    
  res <- signal_tbb[which(signal_tbb$signal==n), ][["interp.res"]]
  time <- index(sim_data[,1])

  tmp <- apply(sim_data, 2, function(x) SpecMTM(MakeEquidistant(time, (coredata(ts(x))-mean(coredata(ts(x)))), dt=res), k=k, nw=nw, detrend=TRUE))
  fit <- lapply(tmp, function(x) SlopeFit(x, 1/max(tscales[[scale]]), 1/min(tscales[[scale]]), bDebug=T))
  scaling[[n]] <- lapply(fit, function(x) list(slope=x$slope, slopesd=x$slopesd))
  names(scaling[[n]]) <- pages.meta$Name  
}

#if(n=="ECHAM5"){
#  next
  #hadCM3 needs to be interpolated or coordinates needs to be obmitted because of NA values
        #len[[i]] <- length(which(is.na(summary$temp[coords[i,]$lon, coords[i,]$lat, ]))) 
        #plot(unlist(len)) maximum gap is 20
        #}
 #   tmp <- apply(sim_data, 2, function(x){
#          #Find the longest consecutive stretch of non-missing values
  #      y <- na.contiguous(zoo(na.approx(as.numeric(coredata(ts(x))), x = time, na.rm = FALSE, maxgap=2), order.by = time))
   #     SpecMTM(MakeEquidistant(index(y), coredata(y)-mean(coredata(y)), dt = res), k=k, nw=nw, detrend=TRUE)
    #    }
     #   )
#  }#
#if(n!="ECHAM5"){
   # }
      