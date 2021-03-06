source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
library(nest)
library(ncdf4)
RMST <- list()

#load data and define parameters: Pages2k
cut_time <- 2020
cut <- FALSE
restrictions <- "loose"
min.res=get_restrictions(restrictions)$min.res
min.range=get_restrictions(restrictions)$min.range
max.hiat= get_restrictions(restrictions)$max.hiat
length.min = get_restrictions(restrictions)$length.min
source("processing/get_rms_pages2k.R")

#hadcrut parameters
hadcrut.len <- 150
source("processing/get_rms_HadCRUT4.R")

#model parameters
meta.res_tbb <- readRDS("helpers/meta_res.Rds")
signal <- signal_tbb %>% filter(type=="model") %>% select(signal)

#compute local mean spectra for specific simulation (i)
i <- seq(1,length(signal$signal),1) #CHOOSE DATA SETS HERE (calculation can require substantial memory and computing power)
for(n in signal[i,]$signal){
    n <- as.character(n)
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
