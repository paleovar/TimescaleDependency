if(!exists("RMST")){
  RMST <- list()
}

#Please download the "PAGES2k_v2.0.0.Rds" from https://springernature.figshare.com/collections/A_global_multiproxy_database_for_temperature_reconstructions_of_the_Common_Era/3285353 into the same directory (`./processing/raw_data`). Our code calls this dataset as `./processing/raw_data/PAGES2k_v2.0.0.Rds`.
pdata <- readRDS("processing/raw_data/PAGES2k_v2.0.0.Rds")

cut_warming_pages <- function(tibble, cut_time, cut, length.min=length.min){
  if(cut==FALSE){return(tibble)}
  if(cut==TRUE){
    res <- tibble %>% mutate(data = purrr::map(data, ~ filter(., time <= cut_time)))
    tmp <- lapply(res$data, function(x) dim(x)[1])
    idx <- which(tmp < length.min)
    if(length(idx)==0){return(res)}
    if(length(idx)!=0){return(res[-which(tmp < length.min),])}
  }
}

#Initialize Meta data
#requires: meta[["Name"]], meta$Lon, meta$Lat, prxlist is a list of zoo objects with index => y BP
prxlist <- list()
meta <- list()
meta.names <- c("ID","Name","Lat","Lon","Elev_masl","Archive",
                "Proxy","Parameter","ParameterInfo","orig.calib.ref",
                "newcalib","Reference","filename","folder","agecol",
                "proxcol","depth","separator","decimal","Calibration",
                "Notes","citekey","NOAA.PANGAEA")
meta.names.ind.pages <- c(1,2,3,4,5,6,7)
ind.pages.temp <- which(sapply(1:length(pdata), function(i) !is.null(pdata[[i]]$paleoData[[1]]$paleoMeasurementTable[[1]]$temperature)))
meta.pages <- data.frame(ID=        sapply(ind.pages.temp, function(i) length(prxlist)+i), #ID
                         Name=      sapply(ind.pages.temp, function(i) names(pdata)[i]), #Name
                         Lat=       sapply(ind.pages.temp, function(i) pdata[[i]]$geo$latitude), #Lat
                         Lon=       sapply(ind.pages.temp, function(i) pdata[[i]]$geo$longitude), #Lon
                         Elev_masl= unlist(sapply(ind.pages.temp, function(i)
                         {
                           y <- pdata[[i]]$geo$elevation
                           if(is.null(y)) return(NA) else return(y)
                         })), #Elev_masl
                         Archive=   sapply(ind.pages.temp, function(i) pdata[[i]]$archiveType), #Archive
                         Proxy=     unlist(sapply(ind.pages.temp, function(i)
                         {
                           y <- pdata[[i]]$paleoData[[1]]$paleoMeasurementTable[[1]]$temperature$proxy
                           if(is.null(y)) return(NA) else return(y)
                         })) #Proxy
)
names(meta.pages) <- meta.names[meta.names.ind.pages]

#convert data to list of zoo objects
prxlist.pages <- lapply(pdata[ind.pages.temp],
                        function(x)
                        {
                          stopifnot(x$paleoData[[1]]$paleoMeasurementTable[[1]]$year$units == "AD");
                          # If your time series has more than one measurement per unique time, 
                          # you're s.o.l and duplicate-time entries will be removed
                          temp <- x$paleoData[[1]]$paleoMeasurementTable[[1]]$temperature$values
                          time <- x$paleoData[[1]]$paleoMeasurementTable[[1]]$year$values
                          count <- table(time)
                          ind.unique <- time %in% as.numeric(rownames(count[count==1]))
                          zoo(temp[ind.unique], time[ind.unique]);
                        })
meta <- rbind(meta[meta.names.ind.pages], meta.pages)
prxlist <- c(prxlist, prxlist.pages)

#exclude proxy that do not match the requirements (see Methods of Ellerhoff & Rehfeld, 2021)
exclude_proxy <- c("Ant-WDC05ABoreholeTreconstruction", "Eur-LakeSilvapla.Larocque-Tobler.2010")
delnames <- exclude_proxy
remid <- c(which(names(prxlist) %in% delnames))

#filter and process the data
if (is.numeric(remid)){
  prxlist<-prxlist[-remid]
  meta<-meta[-remid,]
}

#nest::qualtiy check
check.prxlist <- nest:::quality_check(prxlist, meta, T0=-50, T1=8000, maxhiat=floor(min.res*5), length.min=length.min, min.res=min.res, min.range=min.range)
pages.meta <- check.prxlist$metafilt
prxlist <- prxlist[check.prxlist$ind.ok]
locs <- data.frame(Name=meta$Name[check.prxlist$ind.ok], Lon=pages.meta$Lon, Lat=pages.meta$Lat)

#cut warming trend if necessary
prxtbb <- list_to_tibble(prxlist) %>% rename(Name = model) %>% 
  inner_join(., as_tibble(locs)) %>% inner_join(meta.pages) %>% 
  cut_warming_pages(., cut_time, cut, length.min=length.min)

#filter for hiatuses (i.e. gaps in the data that are too large)
hiat <- lapply(lapply(prxtbb$data, function(x) diff(x$time)), function(x) hiatus=which(abs(x)>=max.hiat))
diffs <- lapply(prxtbb$data, function(x) av.res = mean(diff(x$time)))
names(hiat) <- prxtbb$Name
prxtbb <- prxtbb %>% 
  add_column(interp.res = unlist(diffs))
prxtbb <- prxtbb %>% filter(Name %in% names(hiat[lengths(hiat) == 0L]))
pages.meta <- pages.meta %>% filter(Name %in% names(hiat[lengths(hiat) == 0L]))
prxlist <- prxlist[prxtbb$Name]
tbb <- tibble(Name=character(), data=list())
for(i in seq_along(names(prxlist))){
      tbb <- rbind(tbb, tibble(
      Name = names(prxlist)[[i]],
      data = prxlist[i]
      )
      )
    }
prxlist <- tbb

#compute the spectra
prxtbbspec <- tibble_spec(equidistant(prxtbb), k=3, nw=2)
rm(check.prxlist, locs, meta, meta.pages, pdata, prxlist.pages, hiat, diffs, prxtbb, tbb)
gc()

#compute regional mean spectra for pages2k
lats <- prxtbbspec$Lat
w.lats <- cos(lats*pi/180)/sum(cos(lats*pi/180))
RMST$pages2k <- MeanSpec(tibble_to_list(prxtbbspec), weights=w.lats)
