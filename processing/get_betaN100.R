source("processing/init_processing.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
library(purrr)
save <- T
scale = "cen"
warm_cutoff <- F

#initialization
results <- list()
tscales <- c(200, 10)
min.res = get_min.res(tscales)
min.range = get_min.range(tscales)
max.hiat = get_max.hiat(tscales)
cut_time = 2020

#load data
prxlist <- readRDS("data/proxylist.Rds")
pages.meta <- readRDS("data/pages_meta_scaling.Rds")
specs <- readRDS("data/proxy_spectra.Rds")
scaling <- readRDS("data/scaling_tbb.Rds") %>% filter(model=="pages2k") %>% filter(scale=="cen") %>% rename(Lon=lon, Lat=lat)
df <- prxlist %>% mutate(diffs=purrr::map(data, function(x) diff(index(x)))) %>% select(-data)
tbb <- inner_join(df, specs)
tbb <- inner_join(tbb, scaling)

#Sampling the uncertainties of the scaling coefficients
N=100
results <- replicate(N, 
                 {
    #create artificial proxy with (sample) and without (powerlaw) block averaging
    tbb <- tbb %>% 
      #check if this works
      mutate(tmp = map2(diffs, slope, function(.diffs, .slope) sample_from_proxy(.diffs, (-1)*.slope))) %>% 
      mutate(powerlaw = map2(tmp, diffs, function(.tmp, .diffs) .tmp[["powerlaw"]]), 
             sample = map2(tmp, diffs, function(.tmp, .diffs) .tmp[["sample"]]))
    #sample block averaged proxy to target resolution (ts) and convert powerlaw series to time series (comparison)
    tbb <- tbb %>% mutate(ts = 
                 map2(sample, interp.res, function(.sample, .interp.res) 
                   MakeEquidistant(index(.sample), coredata(.sample), dt=.interp.res))) %>% 
                                    mutate(comparison = 
                                             map2(powerlaw, interp.res, function(.powerlaw, .interp.res)
                                      MakeEquidistant(index(as.ts(.powerlaw)), coredata(as.ts(.powerlaw)), dt = .interp.res)))
    #compute spectrum for block-averaged, equidistant artificial proxy series (ts) and for the comparison 
    tbb <- tbb %>% mutate(spec = 
                                              map2(ts, comparison, function(.ts, .none){SpecMTM(na.remove(.ts), k=k, nw=nw, detrend=TRUE)})) %>% 
                                        mutate(comparison.spec = 
                                                   map2(comparison, ts, function(.comparison, .none){
                                                     SpecMTM(na.remove(.comparison), k=k, nw=nw, detrend=TRUE)}))
    #fit scaling coefficients to spectrum (spec) and comparison spectrum (comparison.spec) and extract standard deviation of fit
    tbb <- tbb %>% mutate(temp = 
                                              map2(spec, comparison.spec, function(.spec, .none) SlopeFit(.spec, 1/max(tscales), 1/min(tscales)))) %>%
                                   mutate(simslope = map2(temp, comparison.spec, function(.temp, .none) .temp$slope)) %>% 
                                  mutate(slopesd = map2(temp, comparison.spec, function(.temp, .none) .temp$slopesd)) %>%
                                  mutate(compslope = 
                                               map2(comparison.spec, temp, function(.comparison.spec, .none) 
                                                 SlopeFit(.comparison.spec, 1/max(tscales), 1/min(tscales))$slope)) %>% 
                                  select(-temp, -tmp)
    tibble(slope=unlist(tbb$simslope), compslopes=unlist(tbb$compslope), slopesd = unlist(tbb$slopesd)) %>% mutate(diffs=abs(compslopes-slope))
    })
  
tmp <- tibble(
    Reduce("+", results[1,])/N, 
    Reduce("+", results[2,])/N,
    Reduce("+", results[3,])/N,
    Reduce("+", results[4,])/N,
    Name = tbb$Name
  ) 
names(tmp) <- c(names(results[,1]), "Name")   
results <- tmp

pages_meta_scaling <- readRDS("data/pages_meta_scaling.Rds")  %>% select(-Elev_masl, -Proxy) %>% rename(long=Lon, lat=Lat)

results <- results %>% left_join(., pages_meta_scaling) %>% add_column(signal="pages2k") %>% select(-compslopes) %>% 
    mutate(slopesd_new = ifelse(signal=="pages2k", sqrt(slopesd**2+diffs**2), slopesd)) %>% 
    select(-diffs) 


bilint_scalingcoeff <- readRDS("processing/bilint_scalingcoeff.Rds") %>% add_column(Archive="model") %>% 
  left_join(., pages_meta_scaling %>% select(-Archive)) %>% mutate(slopesd_new = slopesd) 


beta_N <- rbind(bilint_scalingcoeff , results) %>% add_column(cutoff=cut_time, scale=scale) %>% mutate(Archive=replace(Archive, Archive=="hybrid", "borehole")) %>% 
  mutate(Archive=replace(Archive, Archive=="lake sediment", "lake/marine sediment")) %>% 
  mutate(Archive=replace(Archive, Archive=="marine sediment", "lake/marine sediment")) 

if(save){saveRDS(beta_N, paste0("data/beta_N", N,".Rds"))}
