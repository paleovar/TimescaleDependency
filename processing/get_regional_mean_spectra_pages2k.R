source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
library(nest)

cut_warming_pages <- function(tibble, cut_time=2020, cut=FALSE, length.min=length.min){
  if(cut==FALSE){return(tibble)}
  if(cut==TRUE){
    res <- tibble %>% mutate(data = purrr::map(data, ~ filter(., time <= cut_time)))
    tmp <- lapply(res$data, function(x) dim(x)[1])
    idx <- which(tmp < length.min)
    if(length(idx)==0){return(res)}
    if(length(idx)!=0){return(res[-which(tmp < length.min),])}
  }
}

#load data and define properties
pdata <- readRDS("processing/raw_data/Pages2k_v2.0.0.Rds")
restrictions <- "loose"
min.res=get_restrictions(restrictions)$min.res
min.range=get_restrictions(restrictions)$min.range
max.hiat= get_restrictions(restrictions)$max.hiat
length.min = get_restrictions(restrictions)$length.min

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
n <- length(prxlist)
meta.pages <- data.frame(ID=        sapply(ind.pages.temp, function(i) n+i), #ID
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
pages.prxlist <- prxlist[check.prxlist$ind.ok]
locs <- data.frame(Name=meta$Name[check.prxlist$ind.ok], Lon=pages.meta$Lon, Lat=pages.meta$Lat)

#cut warming trend if necessary
prxtbb <- list_to_tibble(pages.prxlist) %>% rename(Name = model) %>% 
  inner_join(., as_tibble(locs)) %>% inner_join(meta.pages) %>% 
  cut_warming_pages(., length.min=length.min)

#filter for hiatuses (i.e. gaps in the data that are too large)
hiat <- lapply(lapply(prxtbb$data, function(x) diff(x$time)), function(x) hiatus=which(abs(x)>=max.hiat))
diffs <- lapply(prxtbb$data, function(x) av.res = mean(diff(x$time)))
names(hiat) <- prxtbb$Name
prxtbb <- prxtbb %>% 
  add_column(interp.res = unlist(diffs))
prxtbb <- prxtbb %>% filter(Name %in% names(hiat[lengths(hiat) == 0L]))
pages.meta <- pages.meta %>% filter(Name %in% names(hiat[lengths(hiat) == 0L]))

#compute the spectra
prxtbbspec <- tibble_spec(equidistant(prxtbb), k=3, nw=2)
rm(check.prxlist, locs, meta, meta.pages, pdata, prxlist, prxlist.pages, hiat, diffs, prxtbb)
gc()

#compute regional mean spectra for pages2k
RMST <- list()
lats <- prxtbbspec$Lat
w.lats <- cos(lats*pi/180)/sum(cos(lats*pi/180))
RMST$pages2k <- MeanSpec(tibble_to_list(speclist), weights=w.lats)