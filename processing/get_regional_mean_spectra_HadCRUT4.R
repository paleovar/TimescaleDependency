
source("helpers/init.R")
library(ncdf4)

#define parameters
hadcrut.len <- 150

#load data
ncf <- ncdf4::nc_open("processing/raw_data/HadCRUT.4.6.0.0.anomalies.1.nc")

#get variable
#import the data from HadCRUT ncdf files, define global parameters
dat <- list()
dat$lat <- ncdf4::ncvar_get(ncf,varid = "latitude")
dat$lon <- ncdf4::ncvar_get(ncf,varid = "longitude")
dat$tanom <- ncdf4::ncvar_get(ncf,varid = "temperature_anomaly")
dat$time <- ncdf4::ncvar_get(ncf,varid = "time")
dat$time <- as.Date(dat$time,origin="1850-1-1")
dimnames(dat$tanom)<-list(longitude=dat$lon,latitude=dat$lat,time=dat$time)

#get coordinates
coords <- list()
for(i in 1:length(dat$lon)){
  for(j in 1:length(dat$lat)){
    coords[[length(coords)+1]] <- c(dat$lon[i], dat$lat[j])
  }
}

#extract long series
coverage <- function(min.len, dat){
  result <- list()
  for(i in 1:length(coords)){
      loc <- list()
      loc$lon <- coords[[i]][1]
      loc$lat <- coords[[i]][2]
      lat.i <- which.min(abs(dat$lat-loc$lat))
      lon.i <- which.min(abs(dat$lon-loc$lon))
      target.ano <- dat$tanom[lon.i,lat.i,]

      if(all(is.na(target.ano))){
        next
      }

      if(all(is.na(target.ano))==FALSE){
      df <- na.approx(coredata(dat$tanom[lon.i,lat.i,]), x = index(dat$tanom[lon.i,lat.i,]), na.rm = FALSE, maxgap=2)
      y1 <- length(na.contiguous(df))

      if(y1 >= min.len)
        result[[length(result)+1]] <- c(coords[[i]]) #, y1
      }
      }
  return(result)
}

min.length <- 12*hadcrut.len
coords_list <- coverage(min.length, dat)

if(save){
    saveRDS(coords_list, paste0("processing/raw_data/coords/coords_HadCRUT4_len", hadcrut.len, ".Rds"))
}
