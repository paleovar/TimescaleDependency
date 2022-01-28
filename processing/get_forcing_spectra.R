source("processing/init_processing.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
summary <- readRDS("data/forcing_tbb.Rds")
AOD_convert = 20 
k = 3
nw = 2
df.log = 0.02

vol <- summary %>% filter(unit=="AOD") %>% unnest(data) %>% mutate(val =  val*(-AOD_convert), unit="abs dR") %>% group_by(Name, forcing, unit, label) %>% nest() %>% ungroup()

co2 <- summary %>% filter(Name=="co2") %>% unnest(data) %>% mutate(val = 5.35*log(val/278), unit="abs dR") %>% group_by(Name, forcing, unit, label) %>% nest() %>% ungroup()
 
summary <- rbind(summary, vol, co2) %>% unnest(data) %>% rename(dR=val) %>% group_by(Name, forcing, unit, label) %>% nest() %>% ungroup()

spec_sel <- summary %>% mutate(interp.res = as.numeric(purrr::map(data, ~mean(diff(.$time))))) %>% filter(!Name=="ghg_ml") %>%
  filter(unit %in% c("abs dR", "anom dR"), !Name %in% c("lu", "insol")) %>% equidistant() %>% tibble_spec(., k, nw) %>% select(forcing, Name, interp.res, Spec)

spec_sel_list <- spec_sel %>% tibble_to_list() %>% lapply(., function(x) AddConfInterval(x))
names(spec_sel_list) <- spec_sel$Name 

sol_idx <- which(spec_sel$Name %in% c(summary %>% filter(forcing=="sol", Name!="sol") %>% select(Name))$Name)
co2_idx <- which(spec_sel$Name %in% c(summary %>% filter(forcing=="ghg", Name!="ghg") %>% select(Name))$Name)
vol_idx <- which(spec_sel$Name %in% c(summary %>% filter(forcing=="vol", !Name%in% c("CEA", "vol")) %>% select(Name))$Name)

means <- list()
#--------------mean solar forcing----------#
idx <- which(names(spec_sel_list[sol_idx]) == "sol_fro")
selection <- lapply(spec_sel_list[sol_idx][-idx], function(x) cut(x, 500, 2.8))
means$sol <- append(selection, lapply(spec_sel_list[sol_idx][idx], function(x) cut(x, 5.5, 0.001)))
means$sol <- MeanSpec(means$sol)$spec
means$sol <- AddConfInterval(SpecInterpolateSpline(nth_element(means$sol$freq, 1, 3), means$sol))
#LPlot(means$sol)
#--------------mean co2 forcing-----------#
idx <- which(names(spec_sel_list[co2_idx]) == "co2")
selection <- lapply(spec_sel_list[co2_idx][idx], function(x) cut(x, 500, 27)) 
means$co2 <- append(selection, lapply(spec_sel_list[co2_idx][-idx], function(x) cut(x, 100, 0.0013)))
means$co2 <- MeanSpec(means$co2)$spec
freqref <- c(selection$co2$freq, nth_element(means$co2$freq[which.min(abs(means$co2$freq-max(selection$co2$freq))):length(means$co2$freq)], 1, 10))
means$co2 <- SpecInterpolateSpline(nth_element(freqref, 1, 2), means$co2)
means$co2 <- AddConfInterval(lapply(means$co2, function(x) x[-1]))
#LPlot(LogSmooth(means$co2, df.log=2*df.log))
#--------------mean vol forcing-----------#
means$vol <- MeanSpec(spec_sel_list[vol_idx])$spec

spec_sel_list$meansol <- cut(means$sol, 1000, 0.01) 
spec_sel_list$meanco2 <- means$co2
spec_sel_list$meanvol <- means$vol

forc_speclist_smoothed <- lapply(spec_sel_list, function(x) 
LogSmooth(x, df=2*df.log))
forc_speclist_smoothed$meanvol <- cut(forc_speclist_smoothed$meanvol, 500, 0.1)
forc_speclist_smoothed$vol_cro <- cut(forc_speclist_smoothed$vol_cro, 1000, 0.08)
forc_speclist_smoothed$toohey <- cut(forc_speclist_smoothed$toohey, 1000, 0.2)
forc_speclist_smoothed$ghg_ml_hr <- cut(forc_speclist_smoothed$ghg_ml_hr, 1000, 0.0008)
forc_speclist_smoothed$meanco2 <- cut(forc_speclist_smoothed$meanco2, 1500, 0.0015)
forc_speclist_smoothed$meansol <- cut(forc_speclist_smoothed$meansol, 1000, 0.012)

forc_speclist_smoothed_tbb <- list_to_tibble(forc_speclist_smoothed) %>% rename(Name=model) 
if(save){saveRDS(forc_speclist_smoothed_tbb, "data/forcing_spectra.Rds"))}

source("processing/insol.R")
forc_speclist_smoothed <- list()
forc_speclist_smoothed$Berger78 <- result %>% AddConfInterval()
forc_speclist_smoothed_tbb <- list_to_tibble(forc_speclist_smoothed) %>% rename(Name=model) 
if(save){
  saveRDS(forc_speclist_smoothed_tbb %>% unnest(data) %>% filter(freq >= 200) %>% group_by(Name) %>% nest(), "data/forcing_orb_shortterm_spectra.Rds")
  saveRDS(forc_speclist_smoothed_tbb %>% unnest(data) %>% filter(freq <= 10) %>% group_by(Name) %>% nest(), "data/forcing_orb_longterm_spectra.Rds")
}