source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
scaling <- list()
coord <- list()
save <- F #!be careful with overwriting the original files

get_proxies <- F
cut <- F
cut_time <- 2020

scale <-"cen"
if(get_proxies){
    n = "pages2k"
    min.res = get_min.res(tscale[[scale]])
    min.range = get_min.range(tscale[[scale]])
    max.hiat = get_max.hiat(tscale[[scale]])
    length.min = 30
    #filter proxies and compute spectrum
    source("processing/get_rms_pages2k.R")
    if(save){#!be careful with overwriting the original files
      saveRDS(prxlist, "data/proxylist.Rds")
      saveRDS(pages.meta, "data/pages_meta_scaling.Rds")
      saveRDS(prxtbbspec %>% select(-data), "data/proxy_spectra.Rds")
    }
    rm(min.res, min.range, max.hiat, RMST, w.lats)
    specs <- prxtbbspec %>% select(Name, Lon, Lat, Archive, Proxy, interp.res, Spec)
    coords <- data.frame(long=specs$Lon, lat=specs$Lat)
    idx <-c()
    #get scaling coefficients
    for(i in 1:dim(specs)[1]){
      print(i)
      t <- rename(specs$Spec[[i]], spec=temp)
      class(t) <- "spec"
      skip_to_next <- FALSE
      tryCatch(
        fit <- SlopeFit(t, 1/max(tscale[[scale]]), 1/min(tscale[[scale]]))
        , error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) {
        idx <- cbind(idx, i)
        next}  
      
      scaling[[n]][[i]] <- list(slope=fit$slope, slopesd=fit$slopesd)
      
      rm(fit)
    }
    coord[[n]] <- coords
}

#model parameters
meta.res_tbb <- readRDS("helpers/meta_res.Rds")
signal <- signal_tbb %>% filter(type=="model") %>% select(signal) 

#compute local mean spectra for specific simulation (i)
i <- seq(1,length(signal$signal),1) #CHOOSE DATA SETS HERE (calculation can require substantial memory and computing power)
for(n in signal[i,]$signal){
    n <- as.character(n)
    if(n %in% signal_tbb$signal[which(grepl("highres",signal_tbb$signal))]){
      next
    }
    print(n)

    summary <- readRDS(paste0("processing/raw_data/", n, ".Rds"))
    provideDimnames(summary$temp)
    dimnames(summary$temp) <- list(longitude=summary$lon,latitude=summary$lat,time=summary$time)
    coords <- expand.grid(lon=dimnames(summary$temp)$longitude, lat=dimnames(summary$temp)$latitude)
    res <- meta.res_tbb[which(meta.res_tbb$signal==n), ][["interp.res"]]
    
    specs <- list()
    #compute local spectra
    print("specs")
    for(i in 1:dim(coords)[1]){ 
        if(i%%100==0){
          print(i)}
        specs[[i]] <- SpecMTM(MakeEquidistant(summary$time, summary$temp[coords[i,]$lon, coords[i,]$lat, ]-mean(summary$temp[coords[i,]$lon, coords[i,]$lat, ]), dt = res), k=k, nw=nw, detrend=TRUE)
        }
    rm(summary)
    gc()
    #compute scaling coefficients
    print("fits")
    for(i in 1:dim(coords)[1]){
        fit <- SlopeFit(specs[[i]], 1/max(tscale[[scale]]), 1/min(tscale[[scale]]), bDebug=F)
        if(i %in% seq(0, 100000, 100)){
          print(i)
          fit <- SlopeFit(specs[[i]], 1/max(tscale[[scale]]), 1/min(tscale[[scale]]), bDebug=T)
          }
        scaling[[n]][[i]] <- list(slope=fit$slope, slopesd=fit$slopesd)
        rm(fit)
      }
    coord[[n]] <- coords
}

scaling_tbb <- scaling_to_list(scaling, coord, scale)

#once the above for loop was carried out for all i in "signal", the tibble can be saved as
if(save){saveRDS(scaling_tbb, "data/scaling_tbb.Rds")}# !be careful with overwriting the original file
