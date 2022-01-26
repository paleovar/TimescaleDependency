#run get_rms_all.R for all models (listed in signal)
source("processing/get_rms_all.R")
save <- F

RMST <- lapply(RMST, function(x) LogSmooth(x$spec, df=0.02, removeFirst=1, removeLast=floor(3*length(x$freq)/1000)))

for(i in c("MPI-M_highres", "CESM_LM_highres", "ERA5_highres")){
  if(!i %in% names(RMST)){next}
  RMST[[i]] <- cut(RMST[[i]], from=1/2, to=1/(gregorian_year*1.5), index=FALSE)
}

for(n in c("HadCM3", "hadCRUT", "CESM", "ECHAM5")){
  if(!n %in% names(RMST)){next}
  RMST[[n]] <- cut(RMST[[n]], from=2000, to =1/2.5)
}

if(save){
sub <- list()
  for(i in c("CESM_LM", "CESM_LM_highres", "CESM_mean")){
    if(!i %in% names(RMST)){next}
    sub[[i]] <- RMST[[i]]
  }
sub$CESM_LM_mean <- MeanSpec(list(cut(sub$CESM_LM, from=2000, to=1/3, index=FALSE), sub$CESM_LM_highres))$spec
  saveRDS(list_to_tibble(sub), "data/supp/CESM_spectra.Rds")
}

RMST$`MPI-M` <- MeanSpec(list(cut(RMST$`MPI-M`, from=2000, to=1/3, index=FALSE), RMST$`MPI-M_highres`))$spec
RMST$CESM_LM <- MeanSpec(list(cut(RMST$CESM_LM, from=2000, to=1/3, index=FALSE), RMST$CESM_LM_highres))$spec
RMST$ERA5 <- MeanSpec(list(cut(RMST$ERA5, from=2000, to=1/3, index=FALSE), RMST$ERA5_highres))$spec

RMST$CESM_LM_highres <- NULL
RMST$`MPI-M_highres` <- NULL
RMST$ERA5_highres <- NULL

#speclist_smoothed[which(names(speclist_smoothed) %in% c("ERA5", "hadCRUT"))] <- NULL
w <- rep(1, length(names(RMST)))
w[which(names(RMST) %in% names_echam)] <- echam_weights
w[which(names(RMST) %in% control_runs)] <- control_runs_weights

RMST$mean <- MeanSpec(RMST, weights=w)$spec %>% 
  SpecInterpolate(., freqRef = c(RMST$Trace21k$freq, 
                                 RMST$CESM_LM$freq[which.min(abs(RMST$CESM_LM$freq-max(RMST$Trace21k$freq))):length(RMST$CESM_LM$freq)])) %>% 
  SpecInterpolate(., freqRef = nth_element(.$freq, 1, 3)) %>%
  AddConfInterval() %>% cut(., from=10000, to=1/(gregorian_year*1.49))

if(save){
  saveRDS(list_to_tibble(RMST), "data/local_mean_spectra.Rds")) 
}
