source("helpers/init.R")
source("helpers/functions.R")

#---------------------#
#F5
summary <- readRDS("data/forcing_tbb.Rds")
forc_speclist_smoothed_tbb <-  readRDS("data/forcing_spectra_sm.Rds")
cut_colors <- setNames(c("#AA3377", "#228833", "#4477AA", "grey20",  "#CCBB44"), c("sol", "ghg", "vol", "orb", "lu"))

test <- forc_speclist_smoothed_tbb %>% full_join(., summary %>% select(Name, forcing, label))
test$forcing[which(is.na(test$forcing))] <- test$Name[which(is.na(test$forcing))] 

forc_speclist_smoothed_tbb <- test

forc_speclist_smoothed_tbb_orb <- forc_speclist_smoothed_tbb %>% filter(forcing=="Berger78")
forc_speclist_smoothed_tbb <- forc_speclist_smoothed_tbb %>% filter(forcing!="orb") %>% filter(!forcing=="Berger78")

N <- setNames(c(forc_speclist_smoothed_tbb %>% count(forcing) %>% select(n))$n, c(forc_speclist_smoothed_tbb %>% count(forcing) %>% select(forcing))$forcing)

test <- forc_speclist_smoothed_tbb %>% filter(forcing=="meansol") %>% select(data)
idx <- which(forc_speclist_smoothed_tbb$Name == "meansol")
test$data[[1]]$lim.1[[which(test$data[[1]]$lim.1==max(test$data[[1]]$lim.1))]] <- test$data[[1]]$lim.1[[which(test$data[[1]]$lim.1==max(test$data[[1]]$lim.1)) +1]]
test$data[[1]]$lim.2[[which(test$data[[1]]$lim.2==max(test$data[[1]]$lim.2))]] <- test$data[[1]]$lim.2[[which(test$data[[1]]$lim.2==max(test$data[[1]]$lim.2)) +1]]
forc_speclist_smoothed_tbb$data[[idx]]$lim.1 <- test$data[[1]]$lim.1
forc_speclist_smoothed_tbb$data[[idx]]$lim.2 <- test$data[[1]]$lim.2

yrs.period <- rev(c(0.0001, 0.001, 0.01,  0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
yrs.labels <- rev(c(TeX('$10^{-4}$'),TeX('$10^{-3}$'), TeX('$10^{-2}$'), TeX('$10^{-1}$'), TeX('$10^{0}$'), TeX('$10^{1}$'), TeX('$10^{2}$'), TeX('$10^{3}$'), TeX('$10^{4}$'), TeX('$10^{5}$'), TeX('$10^{6}$')))  

plotsol <- ggplot() + 
  theme_td() +
  geom_ribbon(data= forc_speclist_smoothed_tbb %>% filter(forcing=="meansol", Name!="sol") %>%  unnest(data), alpha=0.5, aes(x = 1/freq, ymin=lim.1, ymax=lim.2),  fill="grey60") +
  geom_line(data= forc_speclist_smoothed_tbb %>% filter(forcing=="sol", Name!="sol") %>%  unnest(data), aes(x=1/freq, y = spec, color=label), size=pointsize+0.2) +
  geom_line(data= forc_speclist_smoothed_tbb %>% filter(forcing=="meansol", Name!="sol") %>%  unnest(data), aes(x=1/freq, y = spec), color="black", alpha=0.6,  size=pointsize) +
  scale_color_manual(values=colorRampAlpha(c(cut_colors[["sol"]], "white"), n=N[["sol"]] +1, alpha=1)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)), name=TeX('PSD $S(\\tau)\\, (K^2 yr)$ '), expand=c(0.05, 0.05), limits=c(1e-9, 10), sec.axis = dup_axis(name = NULL, labels = NULL))  +
  scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, limits=c(1e3, 1e-3),  expand=c(0.05, 0.05), name=TeX('period $\\tau\\,(yr)$'), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(legend.position=c(0.7,0.86)) +
  guides(color=guide_legend(ncol=2))
  
print(plotsol)
  
plotghg <- ggplot() +
    theme_td() +
    geom_ribbon(data=forc_speclist_smoothed_tbb %>%filter(forcing=="meanco2")%>%  unnest(data), alpha=0.5, aes(x = 1/freq, ymin=lim.1, ymax=lim.2),  fill="grey60") +
    geom_line(data= forc_speclist_smoothed_tbb %>% filter(forcing=="ghg", Name!="ghg") %>%  unnest(data), aes(x=1/freq, y = spec, color=label), size=pointsize+0.2) +
    geom_line(data= forc_speclist_smoothed_tbb %>%filter(forcing=="meanco2")%>%  unnest(data), aes(x=1/freq, y = spec), color="black", size=pointsize) +
    scale_color_manual(values=colorRampAlpha(c(cut_colors[["ghg"]], "white"), n=N[["ghg"]], alpha=1)) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)), name=TeX('PSD $S(\\tau)\\, (K^2 yr)$ '), limits=c(1e-9, 10), expand=c(0.05, 0.05),  sec.axis = dup_axis(name = NULL, labels = NULL))  +
    scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels,  limits=c(1e3, 1e-3),  expand=c(0.05, 0.05), name=TeX('period $\\tau\\,(yr)$'), sec.axis = dup_axis(name = NULL, labels = NULL)) +
    theme(legend.position=c(0.8,0.8))
  
print(plotghg)
  
plotvol <- ggplot() +
    theme_td() +
    geom_ribbon(data= forc_speclist_smoothed_tbb %>%filter(forcing=="meanvol") %>%  unnest(data), alpha=0.5, aes(x = 1/freq, ymin=lim.1, ymax=lim.2),  fill="grey60") +
    geom_line(data= forc_speclist_smoothed_tbb %>% filter(forcing=="vol", !Name%in% c("CEA", "vol")) %>%  unnest(data), aes(x=1/freq, y = spec, color=label), size=pointsize+0.2) +
    geom_line(data= forc_speclist_smoothed_tbb %>%filter(forcing=="meanvol") %>%  unnest(data), aes(x=1/freq, y = spec), color="black", size=pointsize) +
    scale_color_manual(values=colorRampAlpha(c(cut_colors[["vol"]], "white"), n=N[["vol"]], alpha=1)) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)), name=TeX('PSD $S(\\tau)\\, (K^2 yr)$ '), expand=c(0.05, 0.05), limits=c(1e-9, 10), sec.axis = dup_axis(name = NULL, labels = NULL))  +
    scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, limits=c(1e3, 1e-3), expand=c(0.05, 0.05), name=TeX('period $\\tau\\,(yr)$'), sec.axis = dup_axis(name = NULL, labels = NULL)) +
    theme(legend.position=c(0.8,0.8)) 
  
print(plotvol)

cowplot::plot_grid(
    plotsol  + theme(axis.title.x = element_blank(),
                     legend.position=c(0.48,0.13)), 
    plotghg + theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    legend.position=c(0.7,0.9)), 
    plotvol + theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.title.x = element_blank(),
                    legend.position=c(0.7,0.86)),  
    align = 'h',
    nrow=1,
    label_size = 9,
    labels =c("               (a)", "(b)", "(c)"),
    rel_widths=c(1.25, 1, 1)
  )

#---------------------#
#FS5
cut_colors <- setNames(c("#AA3377", "#228833", "#4477AA", "grey20"), c("meansol", "meanco2", "meanvol", "Berger78"))
forc_speclist_smoothed_tbb  <- bind_rows(forc_speclist_smoothed_tbb, forc_speclist_smoothed_tbb_orb)

plotmeans <- forc_speclist_smoothed_tbb %>% filter(Name %in% c("meansol", "meanvol", "meanco2")) %>% 
  unnest(data) %>%
  ggplot() +
  theme_td() + 
  geom_ribbon(aes(x=1/freq, ymin=lim.1, ymax=lim.2, fill=forcing), alpha=0.3) +
  geom_line(aes(x=1/freq, y = spec, color=forcing), size=pointsize) +
  geom_line(data= forc_speclist_smoothed_tbb_orb %>% unnest(data), aes(x=1/freq, y = spec, color=forcing), size=pointsize) +
  scale_color_manual(values=cut_colors, labels=c("orbital", "CO2", "solar", "volcanic")) +
  guides(fill=FALSE) +
  scale_fill_manual(values=cut_colors) + 
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)), name=TeX('PSD $S(\\tau)\\, (K^2 yr)$ '),  expand=c(0.05, 0.05), limits=c(1e-9, 10), sec.axis = dup_axis(name = NULL, labels = NULL))  +
  scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, limits=c(1e3, 1e-3),  expand=c(0.05, 0.05), name=TeX('period $\\tau\\,(yr)$'), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(legend.position="none") +
  guides(color=guide_legend(ncol=2)) +
  annotate("text", x=c(rep(1000,4)), y=c(5e-7, 5e-8, 5e-6, 5e-9), hjust = 0, label=c("solar", TeX("CO$_2$"), "volcanic", "orbital"), size=9/2.5, color=cut_colors) +
  guides(color = guide_legend(override.aes = list(linetype = 0)))

print(plotmeans)
