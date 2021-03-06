source("helpers/init.R")
source("helpers/functions.R")
#load required packages and initialize data
library(cowplot)
library(forcats)
model <- c(signal_tbb %>% filter(type=="model") %>% select(signal))$signal
obs <- c(signal_tbb %>% filter(type=="obs") %>% select(signal))$signal
signal <- signal_tbb$signal
#load GMST data
GMST <- readRDS("data/GMST_tbb.Rds") %>% filter(!signal %in% c("MPI-M_highres", "CESM_LM_highres", "ERA5_highres", "CESM_LM_cont", "MPI-M_cont"))
GMST <- GMST %>% ungroup() %>% mutate(signal= case_when(signal=="hadCRUT" ~ "hadCRUT4", TRUE ~ signal))

#Compute running mean
GMST <- GMST %>% unnest(data) %>% 
  add_column(year=1*floor(.$time/1)) %>% 
  group_by(signal, year) %>% 
  summarise(median = mean(temp))

#Compute mean from 7 pages2k reconstructions
temp <- GMST %>% filter(signal %in% c("BHM", "CPS", "DA", "M08", "OIE", "PAI", "PCR")) %>% 
  group_by(year) %>% nest() %>% 
  mutate(median = as.numeric(purrr::map(.x=data, .f=  ~mean(.x$median, na.rm = T)))) %>% 
  add_column(signal="pages2k") %>% select(-data) 
GMST <- GMST %>% bind_rows(., temp) %>% filter(!signal %in% c("BHM", "CPS", "DA", "M08", "OIE", "PAI", "PCR"))
GMST <- GMST %>% group_by(signal) %>% mutate(rollm = rollmean(median, 5, align='left', fill=NA))
GMST$median <- coalesce(GMST$rollm, GMST$median)
GMST$rollm <- NULL

#Compute anomalies w.r.t. reference period for HadCM3
temptbb_HadCM3 <- GMST %>% filter(!signal %in% c("ERA5", "pages2k", "hadCRUT4")) %>%  filter(year >= 1800 & year <= 1850) %>% group_by(signal) %>% nest() %>%
  mutate(mean1850=purrr::map(.x = data, .f = ~mean(.x$median, na.rm = T))) %>% select(-data) %>% mutate(mean1850=unlist(mean1850))
modelmean1850 = temptbb_HadCM3 %>% filter(!signal %in% c("HadCM3", "ERA5", "pages2k", "hadCRUT4")) %>% 
  ungroup() %>% select(mean1850) %>% summarise(mean(mean1850)) %>% as.numeric()
temptbb_HadCM3 <- temptbb_HadCM3 %>% filter(signal=="HadCM3")

#Compute anomalies w.r.t. reference period for all other datasets
temptbb <- GMST %>% filter(signal!="HadCM3", year >= 1961 & year <= 1990) %>% group_by(signal) %>% nest() %>%
  mutate(mean6190=purrr::map(.x = data, .f = ~mean(.x$median, na.rm = T))) %>% select(-data) %>% mutate(mean6190=unlist(mean6190))
modelmean6190 = temptbb %>% filter(!signal %in% c("HadCM3", "ERA5", "pages2k", "hadCRUT4")) %>% 
  ungroup() %>% select(mean6190) %>% summarise(mean(mean6190)) %>% as.numeric()
corr <- modelmean6190 - modelmean1850
temptbb_HadCM3$mean6190 <- temptbb_HadCM3$mean1850 + corr
temptbb_HadCM3$mean1850 <- NULL
temptbb <- temptbb %>% rbind(temptbb, temptbb_HadCM3)

#prepare plotting
GMST <- GMST %>% inner_join(., temptbb) %>% inner_join(., signal_tbb) %>% arrange(., type)
levels <- GMST %>% group_by(signal, type) %>% nest() %>% select(signal, type) %>% inner_join(.,signal_tbb) %>% arrange(., order)
cut_colors <-setNames(c(rev(brewer.pal((length(which(levels$signal %in% model))+1), "YlGnBu"))[1:((length(which(levels$signal %in% model))))],  rev(brewer.pal((length(which(levels$signal %in% obs))+1), "OrRd")[2:(length(which(levels$signal %in% obs))+1)])),
                      levels$signal)
#Create main plot
full_plot <- GMST %>% 
  ggplot(aes(y = median-mean6190, x = year, color = signal)) +  
  xlab(TeX('time (yr CE)')) + ylab(TeX('Temperature anomaly (K)')) +
  geom_line(size=pointsize)  + 
  scale_color_manual(values = cut_colors) +
  ylim(c(-1, 1)) +
  theme_td() +
  scale_x_continuous(expand = c(0., 0.), limits=c(0, 2050), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(expand = c(0., 0.), limits=c(-1.1, 1.1), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(legend.position="none")

#get legends
for_model_legend <- GMST %>% 
  filter(signal %in% model) %>%
  ungroup() %>% mutate(signal = factor(signal, levels = unique(.$signal))) %>%
  ggplot(aes(y = median-unlist(mean6190), x = year, color=fct_reorder(signal, order, min))) + 
  geom_line() +
  scale_color_manual(values = cut_colors[names(cut_colors) %in% model], labels= levels$alt_name[levels$type == "model"], name = "Models:") + 
  theme_td() + theme(legend.title = element_text(face="bold"))
for_obs_legend <- GMST %>%
  filter(signal %in% obs) %>%
  ggplot(aes(y = median-unlist(mean6190), x = year, color=signal)) + 
  geom_line() +
  scale_color_manual(values = cut_colors[names(cut_colors) %in% obs], name = "Obs./Recon./Rean.:") +
  theme_td() + theme(legend.title = element_text(face="bold"))

#plot altogether
cowplot::plot_grid(
  full_plot, 
  plot_grid(
    get_legend(for_model_legend), 
    get_legend(for_obs_legend),
    nrow = 2,
    align = "v"
  ), 
  nrow = 1,
  rel_widths = c(2.5,1)
)
