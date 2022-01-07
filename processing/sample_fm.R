if(!exists("samples")){
  samples <- list()
}

summary <- readRDS("data/forcing_tbb.Rds")

samples <- list()
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

samples$volc <- l
rm(l)
gc()

#--------------solar-----------#

spec_sel <- summary %>% filter(forcing=="sol", unit %in% c("abs dR", "anom dR"), Name!= "sol_fro") %>% 
  mutate(interp.res = as.numeric(purrr::map(data, ~mean(diff(.$time))))) %>% 
  equidistant() %>% tibble_spec(., k, nw) %>% select(forcing, Name, interp.res, Spec)

spec_sel_list <- spec_sel %>% tibble_to_list() %>% lapply(., function(x) AddConfInterval(x))
names(spec_sel_list) <- spec_sel$Name 

l <- list()
len <- length(names(spec_sel_list))
for(i in 1:N){
  w <- sample(0:1000, len, replace = TRUE)/1000
  l[[i]] <- MeanSpec(spec_sel_list, weights=w)$spec 
  #l$mean$dof <- NULL
  l[[i]]$lim.1 <- NULL
  l[[i]]$lim.2 <- NULL
  LPlot(l[[i]])
}

samples$sol <- l
rm(l)
gc()


