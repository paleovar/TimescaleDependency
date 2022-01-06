source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
#load data and meta data
GMST <- readRDS("data/GMST_tbb.Rds")
meta.res_tbb <- readRDS("helpers/meta_res.Rds")
save <- F

#process global mean spectra
GMST_equi <- GMST %>% left_join(., meta.res_tbb) %>% mutate(interp.res = replace_na(interp.res, 1)) %>% #cut_warming %>% 
  equidistant %>% select(signal, interp.res, EquiTS)
GMST_spec <- tibble_spec(GMST_equi) %>% select(signal, interp.res, Spec)
speclist <- tibble_to_list(GMST_spec)
names(speclist) <- GMST_spec$signal
pages2k_spectra <- speclist[which(names(speclist) == "BHM"):which(names(speclist) == "PCR")]
speclist[which(names(speclist) == "BHM"):which(names(speclist) == "PCR")] <- NULL
speclist$pages2k <- MeanSpec(pages2k_spectra)$spec
speclist <- lapply(speclist, function(x) PaleoSpec::AddConfInterval(x))

#cut high resolution spectra to merge them with those of lower resolution
for(i in c("MPI-M_highres", "CESM_LM_highres", "ERA5_highres")){
  speclist[[i]] <- cut(speclist[[i]], from=1/2, to=1/(gregorian_year*1.5), index=FALSE)
}

#logsmoothing
speclist_smoothed <- lapply(speclist, function(x) LogSmooth(x, df=0.005, removeFirst=1, removeLast=floor(3*length(x$freq)/1000)))
speclist_smoothed_tbb <- list_to_tibble(speclist_smoothed)
rm(GMST, GMST_equi, GMST_spec, pages2k_spectra, speclist, speclist_tbb)
gc()

#form mean spectra between high and low resolved simulations
speclist_smoothed$`MPI-M` <- MeanSpec(list(cut(speclist_smoothed$`MPI-M`, from=2000, to=1/3, index=FALSE), speclist_smoothed$`MPI-M_highres`))$spec
speclist_smoothed$`MPI-M_highres`  <- NULL

speclist_smoothed$CESM_LM <- MeanSpec(list(cut(speclist_smoothed$CESM_LM, from=2000, to=1/3, index=FALSE), speclist_smoothed$CESM_LM_highres))$spec
speclist_smoothed$CESM_LM_highres  <- NULL

speclist_smoothed$ERA5 <- MeanSpec(list(cut(speclist_smoothed$ERA5, from=2000, to=1/3, index=FALSE), speclist_smoothed$ERA5_highres))$spec
speclist_smoothed$ERA5_highres <- NULL

#exclude some of the spectra from forming the mean (see methods of "Ellerhoff and Rehfeld (2021)")
w <- rep(1, length(names(speclist_smoothed)))
w[which(names(speclist_smoothed) %in% names_echam)] <- echam_weights
w[which(names(speclist_smoothed) %in% control_runs)] <- control_runs_weights

#form mean spectrum
speclist_smoothed$mean <- MeanSpec(speclist_smoothed, weights=w)$spec %>% 
  SpecInterpolate(., freqRef = c(speclist_smoothed$Trace21k$freq, 
                                 speclist_smoothed$CESM_LM$freq[which.min(abs(speclist_smoothed$CESM_LM$freq-max(speclist_smoothed$Trace21k$freq))):length(speclist_smoothed$CESM_LM$freq)])) %>% 
  SpecInterpolate(., freqRef = nth_element(.$freq, 1, 3)) %>%
  AddConfInterval() %>% cut(., from=10000, to=1/(gregorian_year*1.49))
speclist_smoothed_tbb <- list_to_tibble(speclist_smoothed)

#save global mean spectra
if(save){
    saveRDS(speclist_smoothed_tbb, "data/global_mean_spectra.Rds")
}
