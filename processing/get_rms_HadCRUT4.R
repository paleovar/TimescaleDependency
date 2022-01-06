#load data
ncf <- ncdf4::nc_open("processing/raw_data/HadCRUT.4.6.0.0.anomalies.1.nc")

#get variable
#import the data from HadCRUT ncdf files, define global parameters
dat <- list()
dat$lat <- ncdf4::ncvar_get(ncf,varid = "latitude")
dat$lon <- ncdf4::ncvar_get(ncf,varid = "longitude")
dat$tanom <- ncdf4::ncvar_get(ncf,varid = "temperature_anomaly")
gregorian_year <- 365.2425
dat$time <- ncdf4::ncvar_get(ncf,varid = "time")/gregorian_year
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
  result1 <- list()
  result2 <- list()
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
      df <- na.contiguous(df)

      if(y1 >= min.len){
        result1[[length(result1)+1]] <- c(coords[[i]]) #, y1
        result2[[length(result2)+1]] <- MakeEquidistant(t.x=seq(1:length(df))/12, t.y= df-mean(df), dt = 1/12)
        }
      }
      }
  return(list("coords_list"=result1, "temp_list"=result2))
}

seq(1:length(1:10))/12

min.length <- 12*hadcrut.len
result <- coverage(min.length, dat)
length(result$temp_list)

if(save){
    saveRDS(coords_list, paste0("processing/raw_data/coords/coords_HadCRUT4_len", hadcrut.len, ".Rds"))
}

lats <- unlist(lapply(result$coords_list, function(x) x[2]))
w.lats <- cos(lats*pi/180)/sum(cos(lats*pi/180))
specs <- lapply(result$temp_list, function(x) SpecMTM(x, k=3, nw=2, detrend=T))
RMST <- list()
RMST$HadCRUT4 <- MeanSpec(specs, weights=w.lats)
