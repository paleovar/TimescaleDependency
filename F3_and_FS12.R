source("helpers/stacymap.R")
source("helpers/init.R")
pages <- readRDS("data/proxy_scaling.Rds")
tbb <- readRDS("data/scaling_tbb.Rds")
signal_tbb <-readRDS("helpers/signal_tbb.Rds")
model_names <- (signal_tbb %>% filter(type=="model") %>% filter(!signal %in% c("MPI-M_highres", "CESM_LM_highres")) %>% dplyr::select(signal))$signal

plots_cen <- list()
plots <- list()
plots_slopesd <- list()

points_tbb <- pages %>% ungroup() %>% rename(layer=slope) %>% dplyr::select(long, lat, layer, Archive)

points_tbb$layer <- (-1)*points_tbb$layer

points_tbb <- arrange(points_tbb, by=lat)
idx <- which(abs(diff(points_tbb$long)) <5)[which(abs(diff(points_tbb$long)) <5) %in% which(abs(diff(points_tbb$lat))<5)]

points_tbb$long[idx] <- jitter(points_tbb$long[idx], factor=1.5)
points_tbb$lat[idx] <- jitter(points_tbb$lat[idx], factor=1.5)

#####

for(n in model_names){
  subset <- tbb %>% group_by(scale) %>% filter(scale=="cen") %>% filter(model==n) %>% rename(long = lon) 
  
  field_tbb <- subset %>% ungroup() %>% rename(layer=slope) %>% dplyr::select(long, lat, layer)
  field_tbb$layer <- (-1)*field_tbb$layer
  field_tbb <- field_tbb %>% pivot_wider(., names_from = lat, values_from = layer)
  lon <- field_tbb$long
  field_tbb <- field_tbb %>% dplyr::select(-long) %>% as.matrix()
  dimnames(field_tbb) <- list(lon, dimnames(field_tbb)[[2]])
  
  field_sd_tbb <- subset %>% ungroup() %>% rename(layer=slopesd) %>% dplyr::select(long, lat, layer)
  field_sd_tbb <- field_sd_tbb %>% pivot_wider(., names_from = lat, values_from = layer)
  lon1 <- field_sd_tbb$long
  field_sd_tbb <- field_sd_tbb %>% dplyr::select(-long) %>% as.matrix()
  dimnames(field_sd_tbb) <- list(lon1, dimnames(field_sd_tbb)[[2]])
  
  col <- c(rev(brewer.pal(n = 8, name = "Reds")[c(3,8)]), "white", brewer.pal(n = 8, name = "Blues")[c(3)])

  flip <- FALSE
  rot <- FALSE
  if(n %in% c("CESM", "CESM_LM", "CESM_LM_cont", "Trace21k", "Trace21k_orb")){flip <- TRUE}
  if(n %in% c("CESM_LM", "CESM_LM_cont", "MPI-M", "MPI-M_cont", "ECHAM5", "Trace21k", "Trace21k_orb")){rot <- TRUE}
  
  plots_cen[[n]] <- STACYmap(gridlyr = field_tbb, ptlyr = points_tbb, colorscheme =rev(col), ptlyr_shape_var = "Archive",
                graticules = F, box=FALSE, legend_cb = FALSE, legend_num_breaks = 9, flipy = flip, rotate=rot,
                legend_names = list(grid = TeX("$\\beta$"), pt = "Archive")) + 
                ggtitle(as.character(signal_tbb$alt_name[signal_tbb$signal==n])) + 
                theme(plot.title = element_text(size = 9, face="bold", hjust=0.5))
  
  plots[[n]] <- STACYmap(gridlyr = field_tbb, flipy=flip, coastline = FALSE, colorscheme = "temp", rotate=rot)
  plots_slopesd[[n]] <- STACYmap(gridlyr = field_sd_tbb, flipy=flip, coastline = FALSE, colorscheme = "prcp_grd", rotate=rot,
                           graticules = F, box=FALSE, limits=c(0., 0.2)) + 
   ggtitle(as.character(signal_tbb$alt_name[signal_tbb$signal==n])) + 
    theme(plot.title = element_text(size = 9, face="bold", hjust=0.5))
}
#####

prow <- cowplot::plot_grid(
  plots_cen$HadCM3 + theme(legend.position="none", axis.title.x=element_blank(),
                         axis.text.x=element_blank(), plot.margin = unit(c(0.1, 0.1, -0.2, 0.1), "cm")),
  plots_cen$CESM + theme(legend.position="none", axis.title.x=element_blank(),
                            axis.text.x=element_blank(), axis.title.y=element_blank(),
                            axis.text.y=element_blank(), plot.margin = unit(c(0.1, 0.1, -0.2, 0.1), "cm")),
  plots_cen$CESM_LM + theme(legend.position="none", axis.title.x=element_blank(),
                            axis.text.x=element_blank(), axis.title.y=element_blank(),
                            axis.text.y=element_blank(), plot.margin = unit(c(0.1, 0.2, -0.2, 0.1), "cm")),
  plots_cen$CESM_LM_cont + theme(legend.position="none", axis.title.x=element_blank(), plot.margin = unit(c(-0.2, 0.1, 0.2, 0.1), "cm")),
  plots_cen$`MPI-M` + theme(legend.position="none",axis.title.x=element_blank(), axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.1, 0.2, 0.1), "cm"),
                         axis.text.y=element_blank()),
  plots_cen$`MPI-M_cont` + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.2, 0.2, 0.1), "cm"),
                           axis.text.y=element_blank()),
  plots_cen$ECHAM5 + theme(legend.position="none", plot.margin = unit(c(-0.2, 0.1, 0.2, 0.1), "cm")),
  plots_cen$Trace21k_orb + theme(legend.position="none", axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.1, 0.2, 0.1), "cm"),
                         axis.text.y=element_blank()),
  plots_cen$Trace21k + theme(legend.position="none", axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.2, 0.2, 0.1), "cm"),
                           axis.text.y=element_blank()),
  plots_cen$IPSL + theme(legend.position="none", axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.2, 0.2, 0.1), "cm"),
                                 axis.text.y=element_blank()),
  align = 'vh',
  labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)"),
  label_size = 9,
  label_x=0.1,
  label_y=0.98,
  nrow = 4
)

legend1 <- cowplot::get_legend(
  plots_cen$CESM + guides(fill= guide_legend(title=TeX("$\\beta$"), title.vjust = GLOBAL_STACY_OPTIONS$GLOBAL_LEG_TITLE_VJUST,
                                             ncol = 9, direction = 'horizontal', label.position = 'bottom', order=1), 
                          shape=guide_legend(title= TeX("archive"), nrow=1, order=2, title.position="left", vjust=0)) +
                                  theme(legend.box = "vertical",
                              legend.title = element_text(size=11, face="plain"))
)

cowplot::ggdraw() +
  cowplot::draw_plot(prow, 0, 0, 1, 1) +
  cowplot::draw_plot(legend1, 0.42, -0.13, .5, .5)

prow <- cowplot::plot_grid(
  plots_slopesd$HadCM3 + theme(legend.position="none", axis.title.x=element_blank(),
                           axis.text.x=element_blank(), plot.margin = unit(c(0.1, 0.1, -0.2, 0.1), "cm")),   plots_slopesd$CESM + theme(legend.position="none", axis.title.x=element_blank(),
                         axis.text.x=element_blank(), axis.title.y=element_blank(),
                         axis.text.y=element_blank(), plot.margin = unit(c(0.1, 0.1, -0.2, 0.1), "cm")),
  plots_slopesd$CESM_LM + theme(legend.position="none", axis.title.x=element_blank(),
                            axis.text.x=element_blank(), axis.title.y=element_blank(),
                            axis.text.y=element_blank(), plot.margin = unit(c(0.1, 0.2, -0.2, 0.1), "cm")),
  plots_slopesd$CESM_LM_cont + theme(legend.position="none", axis.title.x=element_blank(), plot.margin = unit(c(-0.2, 0.1, 0.2, 0.1), "cm")),
  plots_slopesd$`MPI-M` + theme(legend.position="none",axis.title.x=element_blank(), axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.1, 0.2, 0.1), "cm"),
                            axis.text.y=element_blank()),
  plots_slopesd$`MPI-M_cont` + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.2, 0.2, 0.1), "cm"),
                                 axis.text.y=element_blank()),
  plots_slopesd$ECHAM5 + theme(legend.position="none", plot.margin = unit(c(-0.2, 0.1, 0.2, 0.1), "cm")),
  plots_slopesd$Trace21k_orb + theme(legend.position="none", axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.1, 0.2, 0.1), "cm"),
                                 axis.text.y=element_blank()),
  plots_slopesd$Trace21k + theme(legend.position="none", axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.2, 0.2, 0.1), "cm"),
                             axis.text.y=element_blank()),
  plots_slopesd$IPSL + theme(legend.position="none", axis.title.y=element_blank(),plot.margin = unit(c(-0.2, 0.2, 0.2, 0.1), "cm"),
                         axis.text.y=element_blank()),
  align = 'vh',
  labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)"),
  label_size = 9,
  label_x=0.1,
  label_y=0.98,
  nrow = 4
)

legend1 <- cowplot::get_legend(
  plots_slopesd$CESM + guides(fill= guide_legend(title=TeX("std. $\\Delta \\beta$"),
                                             direction = 'horizontal', label.position = 'bottom', order=1)) +
    theme(legend.box = "vertical",
          legend.title = element_text(size=11, face="plain"))
)

cowplot::ggdraw() +
  cowplot::draw_plot(prow, 0, 0, 1, 1) +
  cowplot::draw_plot(legend1, 0.42, -0.13, .5, .5)

