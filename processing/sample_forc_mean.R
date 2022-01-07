#--------------mean vol forcing-----------#
rm(list = ls(all.names = TRUE)) 
gc()

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
dir <- sub('paper*', '', dir)

source(paste0(dir, "/scripts/0_init.R"))
source(paste0(dir, "/scripts/0_helpfunctions.R"))
summary <- readRDS(paste0(dir, "data/forcings/summary.Rds"))

l <- list()
len <- nrow(summary %>% filter(unit=="AOD"))

N <- 100

for(i in 1:N){
  aod_convert <- sample(-18:-25, 1, replace = TRUE)
  w <- sample(0:1000, len, replace = TRUE)/1000
  vol <- summary %>% filter(unit=="AOD") %>% unnest(data) %>% mutate(val =  val*(-aod_convert), unit="abs dR") %>% group_by(Name, forcing, unit, label) %>% nest() %>% ungroup()
  vol <- vol %>% unnest(data) %>% rename(dR=val) %>% group_by(Name, forcing, unit, label) %>% nest() %>% ungroup()
  spec_sel <- vol %>% mutate(interp.res = as.numeric(purrr::map(data, ~mean(diff(.$time))))) %>% filter(!Name=="ghg_ml") %>%
    filter(unit %in% c("abs dR", "anom dR"), !Name %in% c("lu", "insol")) %>% equidistant() %>% tibble_spec(., k, nw) %>% select(forcing, Name, interp.res, Spec)
  spec_sel_list <- spec_sel %>% tibble_to_list() %>% lapply(., function(x) AddConfInterval(x))
  names(spec_sel_list) <- spec_sel$Name 
  spec_sel_list <- lapply(spec_sel_list, function(x) cut(x, from=5000, 12/9,index=FALSE))
  l[[i]] <- MeanSpec(spec_sel_list, weights=w)$spec 
    #l$mean$dof <- NULL
  l[[i]]$lim.1 <- NULL
  l[[i]]$lim.2 <- NULL
  LPlot(l[[i]])
}

if(save){
saveRDS(l, paste0(dir, "/paper/results/sampling/sample_vol_spec_N", N, ".Rds"))
}
#--------------co2+solar-----------#

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
dir <- sub('paper*', '', dir)

source(paste0(dir, "/scripts/0_init.R"))
source(paste0(dir, "/scripts/0_helpfunctions.R"))
summary <- readRDS(paste0(dir, "data/forcings/summary.Rds"))

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

#--------------mean solar forcing-----------#
idx <- which(names(spec_sel_list[sol_idx]) == "sol_fro")
selection <- spec_sel_list[sol_idx][-idx]
#selection <- append(selection, lapply(spec_sel_list[sol_idx][idx], function(x) cut(x, 5.5, 0.001)))
#means$sol <- selection
#selection$sol_fro

sel_tbb <- list_to_tibble(selection)

ggplot() +
  theme_physrev(textsize) +
  geom_ribbon(data= sel_tbb %>%  unnest(data), alpha=0.5, aes(x = 1/freq, ymin=lim.1, ymax=lim.2,  fill=model)) +
  geom_line(data= sel_tbb %>%  unnest(data), aes(x=1/freq, y = spec, color=model), size=pointsize+0.2) +
  # geom_line(data= forc_speclist_smoothed_tbb %>% filter(forcing=="meansol", Name!="sol") %>%  unnest(data), aes(x=1/freq, y = spec), color="black", alpha=0.6,  size=pointsize) +
  #  scale_color_manual(values=colorRampAlpha(c(cut_colors[["sol"]], "white"), n=N[["sol"]] +1, alpha=1)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)), name=TeX('power spectral density ($W^2 m^{-4} yr$)'), expand=c(0.05, 0.05), limits=c(1e-9, 10), sec.axis = dup_axis(name = NULL, labels = NULL))  +
  scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, limits=c(1e3, 1e-3),  expand=c(0.05, 0.05), name=TeX('time period ($yr$)'), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(legend.position=c(0.7,0.86)) #+
# guides(color=guide_legend(ncol=2))

l <- list()
len <- length(names(selection))
for(i in 1:N){
  w <- sample(0:1000, len, replace = TRUE)/1000
  l[[i]] <- MeanSpec(selection, weights=w)$spec 
  #l$mean$dof <- NULL
  l[[i]]$lim.1 <- NULL
  l[[i]]$lim.2 <- NULL
  LPlot(l[[i]])
}

if(save){
saveRDS(l, paste0(dir, "/paper/results/sampling/sample_solar_spec_N", N, ".Rds"))
}
