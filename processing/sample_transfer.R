source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")

#be patient with the following lines: sampling of spectra
samples <- list()
source("processing/sample_fm.R")
source("processing/sample_rm.R")
source("processing/sample_mm.R")

#load and process data
forc_co2_orb <- rbind(
    readRDS("data/forcing_spectra.Rds") %>% filter(Name == "meanco2") %>% 
    unnest(data) %>% filter(freq < 9/12) %>% group_by(Name) %>% nest() %>% ungroup()
    , readRDS("data/forcing_orb_longterm_spectra.Rds") %>% filter(Name == "Berger78") %>% 
    unnest(data) %>% filter(freq < 9/12) %>% group_by(Name) %>% nest() %>% ungroup() %>% mutate(Name = ifelse(Name == "Berger78", "meanorb", Name))
    )

comp_forc <- lapply(as.list(forc_co2_orb)$data, function(x){
  x$forcing <- NULL
  x$label <- NULL
  class(x) <- "spec"
  return(x)
}
)
names(comp_forc) <- forc_co2_orb$Name

forc_sol <- samples$sol
forc_vol <- samples$vol

resp_recons <- samples$recons

resp_recons2 <- list()
resp_recons2$ERA5 <- tibble_to_list(readRDS("data/global_mean_spectra.Rds") %>% filter(signal=="ERA5"), name.data="data")[[1]]
resp_recons2$ERA5 <- cut(resp_recons2$ERA5, from=5, to=length(resp_recons2$ERA5$freq), index=TRUE)
resp_recons2$hadCRUT <-  tibble_to_list(readRDS("data/global_mean_spectra.Rds") %>% filter(signal=="hadCRUT4"), name.data="data")[[1]]
resp_recons2$hadCRUT <- cut(resp_recons2$hadCRUT, from=13, to=length(resp_recons2$hadCRUT$freq), index=TRUE)

resp_models_wENSO <- samples$models_wENSO
resp_models_woENSO <- samples$models_woENSO

#Sampling of the mean spectra and the gain 
N_sample <- 50 #be carful with changing this parameter. It requires a lot of computing time to go to a larger number of samples
sample_raw <- list()
#GAIN
for(transfertarget in c("recons","forc", "models_woENSO", "models_wENSO")){
    l <- list()
    print(transfertarget)
    for(i in 1:N_sample){
        specList <- list()
        if(transfertarget == "forc"){
        #forcing
        comp_forc$sol <- forc_sol[[sample(1:length(forc_sol), 1, replace = TRUE)]]
        comp_forc$vol <- forc_vol[[sample(1:length(forc_vol), 1, replace = TRUE)]]
        
        specList$forc <- MeanSpectrum(comp_forc)$spec
        
        specList$forc$spec <- specList$forc$spec * length(names(comp_forc))
        specList$forc$lim.1 <-  specList$forc$lim.1 * length(names(comp_forc))
        specList$forc$lim.2 <- specList$forc$lim.2 * length(names(comp_forc))
        }
        
        #GMST
        if(transfertarget == "recons"){
        idx1 <- sample(1:length(resp_recons), 1, replace=TRUE)
        resp_recons2$pages <- resp_recons[[idx1]][[sample(1:length(resp_recons[[idx1]]), 1, replace=TRUE)]]
        specList$resp <- MeanSpec(resp_recons2)$spec
        }
        
        if(transfertarget == "models_woENSO"){
        specList$resp <- resp_models_woENSO[[sample(1:length(resp_models_woENSO), 1, replace=TRUE)]]
        }
        if(transfertarget == "models_wENSO"){
        specList$resp <- resp_models_wENSO[[sample(1:length(resp_models_wENSO), 1, replace=TRUE)]]
        }
        
        #transfer
        l[[i]] <- specList[[1]]
        l[[i]]$lim.1 <- NULL
        l[[i]]$lim.2 <- NULL
        #LPlot(l[[i]])
        
        comp_forc$sol <- NULL
        comp_forc$vol <- NULL
        idx1 <- NULL
    }

    l <- lapply(l, function(x) cut(x, 500, 2, index=FALSE))
    names(l) <- as.character(seq(1:length(t)))

    sample_raw[[transfertarget]] <- list_to_tibble(l) %>% unnest(data) %>% group_by(freq)  %>% 
    summarize(
    m = mean(spec, na.rm=TRUE),  
    sd_up = quantile(spec, c(0.95), na.rm=TRUE),
    sd_down = quantile(spec, c(0.05), na.rm=TRUE)
    )
}

N_sample <- 50 #be carful with changing this parameter. It requires a lot of computing time to go to a larger number of samples
sample_transfer <- list()
#SPECTRA
for(transfertarget in c("recons","forc", "models_woENSO", "models_wENSO")){
    l <- list()
    print(transfertarget)
    for(i in 1:N_sample){
        if(transfertarget=="forc"){next}
    specList <- list()
    #forcing
    comp_forc$sol <- forc_sol[[sample(1:length(forc_sol), 1, replace = TRUE)]]
    comp_forc$vol <- forc_vol[[sample(1:length(forc_vol), 1, replace = TRUE)]]
    
    specList$forc <- MeanSpectrum(comp_forc)$spec
    
    specList$forc$spec <- specList$forc$spec * length(names(comp_forc))
    specList$forc$lim.1 <-  specList$forc$lim.1 * length(names(comp_forc))
    specList$forc$lim.2 <- specList$forc$lim.2 * length(names(comp_forc))

    #GMST
    if(transfertarget == "recons"){
        idx1 <- sample(1:length(resp_recons), 1, replace=TRUE)
        resp_recons2$pages <- resp_recons[[idx1]][[sample(1:length(resp_recons[[idx1]]), 1, replace=TRUE)]]
        specList$resp <- MeanSpec(resp_recons2)$spec
        }
    
    if(transfertarget == "models_woENSO"){
        specList$resp <- resp_models_woENSO[[sample(1:length(resp_models_woENSO), 1, replace=TRUE)]]
    }
    if(transfertarget == "models_wENSO"){
        specList$resp <- resp_models_wENSO[[sample(1:length(resp_models_wENSO), 1, replace=TRUE)]]
    }
    
    #transfer
    l[[i]] <- transferSpec(specList, input=1, output=2)$spec
    l[[i]]$lim.1 <- NULL
    l[[i]]$lim.2 <- NULL
    #LPlot(l[[i]])
    
    comp_forc$sol <- NULL
    comp_forc$vol <- NULL
    idx1 <- NULL
    }

    l <- lapply(l, function(x) cut(x, 500, 2, index=FALSE))
    names(l) <- as.character(seq(1:length(t)))

    sample_transfer[[transfertarget]] <- list_to_tibble(l) %>% unnest(data) %>% group_by(freq)  %>% summarize(
    m = mean(spec, na.rm=TRUE),  
    sd_up = quantile(spec, c(0.95), na.rm=TRUE),
    sd_down = quantile(spec, c(0.05), na.rm=TRUE)
    )
}

for(i in names(sample_raw)){
    sample_raw[[i]] <- sample_raw[[i]] %>% add_column(name=i)
}
for(i in names(sample_transfer)){
    sample_transfer[[i]] <- sample_transfer[[i]] %>% add_column(name=i)
}

#rbind

if(save){
    saveRDS(sample_raw, paste0("data/transfer_", N_sample, ".Rds"))
    saveRDS(sample_transfer, paste0("data/sample_spec_", N_sample, ".Rds"))
}

