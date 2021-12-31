
source("helpers/init.R")
source("helpers/functions.R")
meta.res_tbb <- readRDS("helpers/meta_res.Rds")
save <- F

signal <- signal_tbb %>% filter(type=="model") %>% select(signal)
k=3
nw=2

RMST <-list()
for(n in signal[2,]){
    speclist <- list()
    speclist_lat <- list()    
    print(n)
    if(n %in% c("hadCRUT4", "pages2k")){next} #will be treated separately

    summary <- readRDS(paste0("processing/raw_data/", n, ".Rds"))
    provideDimnames(summary$temp)

    dimnames(summary$temp) <- list(longitude=summary$lon,latitude=summary$lat,time=summary$time)
    coords <- expand.grid(lon=dimnames(summary$temp)$longitude, lat=dimnames(summary$temp)$latitude)
    res <- meta.res_tbb[which(meta.res_tbb$signal==n), ][["interp.res"]]
    speclist <- list()
    for(i in 1:dim(coords)[1]){
        if(i%%100==0){print(i)}
        speclist[[i]] <- SpecMTM(MakeEquidistant(summary$time, summary$temp[coords[i,]$lon, coords[i,]$lat, ]-mean(summary$temp[coords[i,]$lon, coords[i,]$lat, ]), dt = res), k=k, nw=nw, detrend=TRUE)
        }
    rm(summary)
    gc()

    #latitudinal mean spectra
    coords$index <- index(coords)
    for (i in levels(coords$lat)){
      if (i %in% coords$lat){
      idx <- coords[which(coords$lat == i),]$index
      speclist_lat[[i]] <- MeanSpec(speclist[idx])
      }
    }
  
    #area-weighted local mean spectra
    lats <- as.numeric(names(speclist_lat))
    w.lats <- cos(lats*pi/180)/sum(cos(lats*pi/180))
    speclist_lat <- lapply(speclist_lat, function(x) x$spec)
    RMST[[n]] <- MeanSpec(speclist_lat, weights=w.lats)
    rm(speclist_lat, speclist)
    gc()
}
