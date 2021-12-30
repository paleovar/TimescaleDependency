source("helpers/init.R")
source("helpers/functions.R")
#load required packages and data
library(cowplot)
library(grid)
library(gridExtra)
summary <- readRDS("data/forcing_tbb.Rds") 

#define plotting parameters
cut_colors <- setNames(c("#AA3377", "#228833", "#4477AA", "grey20"), c("sol", "ghg", "vol", "orb"))
N <- setNames(c(summary %>% count(forcing) %>% select(n))$n, c(summary %>% count(forcing) %>% select(forcing))$forcing)
forc <- unique(summary$forcing)
ylab_forc_s <- setNames(c(TeX("RF $(W/m^{2})$"), TeX(" conc. $(ppm)$"), TeX("AOD"), TeX(" RF $(W/m^{2})$")), c("sol", "ghg", "vol", "orb"))
leg_forc <- setNames(c(TeX("TSI"), TeX("CO$_2$"), TeX("volcanic"), TeX("insolation $65°$N")), c("sol", "ghg", "vol", "orb"))
pointsize=0.3
plots <- list()

#Compute means
summary_yearly <- summary %>%  filter(!Name %in% c("sol_fro", "sol", "vol_cro", "ghg", "n2o", "ch4", "vol", "ghg_ml", "ghg_ml_hr"))
tmp <- summary_yearly %>% mutate(interp.res = as.numeric(purrr::map(data, ~mean(diff(.$time)))))
for(i in which(tmp$interp.res != 1)){
  summary_yearly$data[[i]] <- tmp$data[[i]] %>%  add_column(year=1*floor(.$time/1)) %>% group_by(year) %>% 
  summarise(median = mean(val)) %>% rename(val=median, time=year)
}

#create subplots for each forcing type
if(all(forc %in% names(cut_colors))){
for(n in forc){
  t <- theme(panel.background = element_rect(fill = "transparent",colour = NA),
             plot.background = element_rect(fill = "transparent",colour = NA),
             legend.title = element_blank(),
             legend.background = element_blank(),
             legend.box.background = element_blank(),
             legend.text = element_text( size=9),
             axis.ticks.y = element_line(colour=cut_colors[[n]]),
             axis.ticks.length = unit(-1.4, "mm"),
             axis.title.y =  element_text(colour=cut_colors[[n]]), 
             axis.title.x=element_blank(), 
             legend.key.size = unit(0.2, "cm")) 
  if(n =="ghg"){
  t <-  t +  theme(legend.position = c(0.18,0.3),
                 axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm"), colour=cut_colors[[n]]),
                 axis.text.x = element_blank(),
                 axis.line.x=element_blank(),
                 axis.ticks.x=element_blank()) 
  }
  if(n =="sol"){
    t <-  t +  theme(legend.position = c(0.44,-0.1),
                     axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm"), colour=cut_colors[[n]]),
                     axis.text.x = element_blank(),
                     axis.line.x=element_blank(),
                     axis.ticks.x=element_blank()) 
  }
  if(n == "vol"){
    t <- t + theme(legend.position = c(0.18,0.6),
                axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm"), colour=cut_colors[[n]]),
                axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"), colour="black"),
                axis.line.x.bottom = element_blank(), 
                axis.ticks.x.bottom = element_blank(), 
                axis.text.x.bottom = element_blank()) 
  }
  if(n == "orb"){
    t <-   t + theme(legend.position = c(0.17,0.25),
                axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm"), colour=cut_colors[[n]]),
                axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"), colour="black"),
                axis.line.x.top = element_blank(), 
                axis.ticks.x.top = element_blank(), 
                axis.text.x.top= element_blank()) 
  }
  
  plots[[n]] <- summary_yearly %>% filter(forcing==n) %>%  unnest(data)%>%
    ggplot(aes(x = time, y = val, color=label)) +
    geom_line(size=pointsize) + scale_color_manual(values=colorRampAlpha(c(cut_colors[[n]], "white"), n=N[[n]] +1, alpha=1)) +
    theme_td() + labs(y=ylab_forc_s[[n]]) +
    theme_classic(base_size=9) + 
    scale_x_continuous(expand=c(0.,0), sec.axis = dup_axis(name = NULL, labels = NULL), limits=c(800, 2020)) +
    scale_y_continuous(expand=c(0.,0.), sec.axis = dup_axis(name = NULL, labels = NULL)) +
    t 
}
}

#Create full plot
plt <- cowplot::plot_grid(plots$vol + 
                    scale_y_continuous(expand=c(0.,0.), sec.axis = dup_axis(name = NULL, labels = NULL), limits=c(0, 1.1)) + 
                    scale_color_manual(values=colorRampAlpha(c(cut_colors[["vol"]],"white"), n=N[["vol"]], alpha=1)), 
                  plots$sol + scale_color_manual(values=colorRampAlpha(c(cut_colors[["sol"]], "white"), n=N[["sol"]], alpha=1)) + 
                    guides(color=guide_legend(ncol=3)), 
                  plots$ghg,
                  plots$orb + scale_y_continuous(expand=c(0.,0.), sec.axis = dup_axis(name = NULL, labels = NULL), limits=c(478.6, 481.75)), 
                  align="v", 
                  ncol=1)
y.top <- grid::textGrob(TeX(' '), 
                            gp=gpar(col="black", fontsize=9))
x.grob <- grid::textGrob(TeX('time $(yr\\,CE)$'), 
                   gp=gpar(col="black", fontsize=9))
#add labels plot
plotgrd <- gridExtra::grid.arrange(arrangeGrob(plt, bottom=x.grob, top=y.top, right=y.top))
#plot altogether
ggdraw(plotgrd) + 
  draw_label(TeX("solar"), x = 0.98, y = 0.63, color=cut_colors[["sol"]], angle=-90,  size=9) +
  draw_label(TeX("volcanic"), x = 0.98, y = 0.88, color=cut_colors[["vol"]], angle=-90,  size=9) +
  draw_label(TeX("CO$_2$"), x = 0.98, y = 0.4, color=cut_colors[["ghg"]], angle=-90,  size=9)  +
  draw_label(TeX("insolation $65°$N"), x = 0.98, y = 0.2, color=cut_colors[["orb"]], angle=-90,  size=9) 

#repeat the same for highly resolved radiative forcing
summary_highres <- summary %>% filter(Name %in% c("vol_cro", "sol_fro", "ghg_ml"))
add_forc <- unique(summary_highres$forcing)
plots <- list()
for(n in add_forc){
  if(n=="vol"){next}
    t <-  theme(legend.position = c(0.85,0.2),
                panel.background = element_rect(fill = "transparent",colour = NA),
                plot.background = element_rect(fill = "transparent",colour = NA),
                legend.title = element_blank(),
                legend.background = element_blank(),
                legend.box.background = element_blank(),
                legend.text = element_text(color=cut_colors[[n]], size=9),
                axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm"), colour=cut_colors[[n]]),
                axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"), colour="black"),
                axis.ticks.y = element_line(colour=cut_colors[[n]]),
                axis.ticks.length = unit(-1.4, "mm"),
                axis.title.x=element_blank(),
                axis.title.y=element_text(colour=cut_colors[[n]]))
  if(n=="sol"){t <- t + theme(axis.text.x.bottom=element_blank(), axis.ticks.x.bottom = element_blank(), axis.line.x.bottom=element_blank())  }
  if(n=="ghg"){t <- t + theme(axis.text.x.top=element_blank(), axis.ticks.x.top = element_blank(), axis.line.x.top=element_blank())  }
    
  plots[[n]] <- summary_highres %>% filter(forcing==n) %>% unnest(data)%>%
    ggplot(aes(x = time, y = val, color=label)) +
    geom_line(size=pointsize) + scale_color_manual(values=colorRampAlpha(c(cut_colors[[n]], "white"), n=1, alpha=1)) +
    theme_td() + labs(y=ylab_forc_s[[n]]) +
    theme_classic(base_size=9) + 
    scale_x_continuous(expand=c(0.,0), sec.axis = dup_axis(name = NULL, labels = NULL), limits=c(1970, 2020)) +
    scale_y_continuous(expand=c(0.,0.), sec.axis = dup_axis(name = NULL, labels = NULL)) +
    t +   guides(color = guide_legend(override.aes = list(linetype = 0)))
}
y.grob <- textGrob(TeX('radiative forcing $(W/m^{2})$'), 
                   gp=gpar(col="black", fontsize=9), rot=90)
y.top <- textGrob(TeX(' '), 
                            gp=gpar(col="black", fontsize=9))
x.grob <- textGrob(TeX('time $(yr\\,CE)$'), 
                   gp=gpar(col="black", fontsize=9))
plt <- plot_grid(plots$sol + scale_y_continuous(expand=c(0.,0.), sec.axis = dup_axis(name = NULL, labels = NULL), limits=c(1362, 1368.5)), plots$ghg, align="v", ncol=1)
plotgrd <- grid.arrange(arrangeGrob(plt, bottom=x.grob, top=y.top, right=y.top))
#plot altogether
ggdraw(plotgrd) + 
  draw_label(TeX("solar"), x = 0.98, y = 0.75, color=cut_colors[["sol"]], angle=-90,  size=9) +
  draw_label(TeX("CO$_2$"), x = 0.98, y = 0.31, color=cut_colors[["ghg"]], angle=-90,  size=9)
