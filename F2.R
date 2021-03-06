source("helpers/init.R")
source("helpers/functions.R")
library(cowplot)
library(forcats)

#load spectra
gmspectra <- readRDS("data/global_mean_spectra.Rds")
lmspectra <- readRDS("data/local_mean_spectra.Rds")

model <- c(signal_tbb %>% filter(type=="model"))$signal
obs <- c(signal_tbb %>% filter(type=="obs"))$signal

#get mean spectrum and filter spectra for global mean surface temperature
mean_df <- gmspectra %>% filter(signal == "mean") %>% unnest(data)
gmspectra <- gmspectra %>% filter(!signal %in% c("MPI-M_cont", "CESM_LM_cont"))  %>%
  inner_join(., signal_tbb) %>% 
  arrange(., signal) %>% arrange(., type)

#assign colors for plotting
levels <- gmspectra %>% group_by(signal, type)  %>% select(-data)  %>% arrange(., order)
cntlevels <- levels %>% group_by(type) %>% count()
cut_colors <-setNames(c(rev(brewer.pal(cntlevels[which(cntlevels$type=="model"),]$n+1, "YlGnBu"))[1:cntlevels[which(cntlevels$type=="model"),]$n], rev(brewer.pal(cntlevels[which(cntlevels$type=="obs"),]$n +1, "OrRd")[2:(cntlevels[which(cntlevels$type=="obs"),]$n+1)]), "black"), c(levels$signal, "mean"))

#create plot for global mean spectra
full_plot1 <- plot_spec(gmspectra, ylims=c(1e-8, 2e3), xlims=c(1e3, 1e-3)) +
  geom_ribbon(data=mean_df, alpha=0.1, aes(x = 1/freq, ymin=lim.1, ymax=lim.2, fill="mean"), linetype=0) +
  geom_line(data=mean_df, aes(x=1/freq,y=spec), color='black',size=0.4) +
  geom_abline(intercept = c(-5, -0.5), slope = c(-1, -2), lty=2, size=0.2) +
  guides(fill=FALSE) + scale_fill_manual(values=cut_colors) + 
  scale_color_manual(values = cut_colors) + 
  theme(legend.position="none") +
  annotate("text", x = c(0.02, 20), y=c(0.01, 0.00003), label = c("β=2", "β=1"), size=notationsize) +
  annotate("text", x = c(350), y=c(1e-7), label = c("global mean"), color=("black"), size=notationsize) +
  annotate("text",x = c(500), y=c(500), label = c("(b)"), size=notationsize)

#get legend for model simulations
for_model_legend1 <- gmspectra %>% unnest(data) %>%
  filter(signal %in% model) %>%
  ungroup() %>% mutate(signal = factor(signal, levels = unique(.$signal))) %>%
  ggplot(aes(x = 1/freq, y = spec, ymin=lim.1, ymax=lim.2,  fill=fct_reorder(signal, order, min))) + 
  geom_line(aes(color=fct_reorder(signal, order, min))) +
  scale_color_manual(values = cut_colors, labels= levels$alt_name[levels$type == "model"], name = "model simulations:") + 
  theme(legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.65, 'lines'),
        legend.title = element_text(size=9*0.75, face="bold"),
        legend.text = element_text(size=9*0.75))

#get legend for observational data
for_obs_legend1 <- gmspectra %>% unnest(data) %>%
  filter(signal %in% obs) %>%
  ggplot(aes(x = 1/freq, y = spec, ymin=lim.1, ymax=lim.2,  fill=signal)) + 
  geom_line(aes(color=signal)) +
  scale_color_manual(values = cut_colors, name = "observation-based:") + 
  theme(legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.65, 'lines'),
        legend.title = element_text(size=9*0.75, face="bold"),
        legend.text = element_text(size=9*0.75))

#get local mean spectra and filtering
lmspectra_mean_df <- lmspectra %>% filter(signal == "mean") %>% unnest(data)
lmspectra <- lmspectra %>% filter(!signal %in% c("MPI-M_cont", "CESM_LM_cont")) %>% inner_join(., signal_tbb) %>% arrange(., signal) %>% arrange(., type)

#assign colors for plotting
levels <- lmspectra %>% group_by(signal, type) %>% select(-data) %>% arrange(., order)
cntlevels <- levels %>% group_by(type) %>% count()
cut_colors <-setNames(c(rev(brewer.pal(cntlevels[which(cntlevels$type=="model"),]$n+1, "YlGnBu"))[1:cntlevels[which(cntlevels$type=="model"),]$n], rev(brewer.pal(cntlevels[which(cntlevels$type=="obs"),]$n +1, "OrRd")[2:(cntlevels[which(cntlevels$type=="obs"),]$n+1)]), "black"), c(levels$signal, "mean"))

#create plot for local mean spectra
full_plot2 <- plot_spec(lmspectra, ylims=c(1e-5, 2e3), xlims=c(1e3, 1e-3)) +
  geom_ribbon(data=lmspectra_mean_df, alpha=0.3, aes(x = 1/freq, ymin=lim.1, ymax=lim.2, fill="mean"), linetype=0) +
  geom_abline(intercept = c(-3, 2.5), slope = c(-1, -2), lty=2, size=0.2) +
  geom_line(data=lmspectra_mean_df, aes(x=1/freq,y=spec), color='black',size=0.4) +
  guides(fill=FALSE) + scale_fill_manual(values=cut_colors) + 
  scale_color_manual(values = cut_colors) + 
  theme(legend.position="none") +
  annotate("text", x = c(0.02, 5), y=c(10, 0.0003), label = c("β=2", "β=1"), size=notationsize) +
  annotate("text",x = c(500), y=c(500), label = c("(a)"), size=notationsize) +
  annotate("text",x = c(350), y=c(1e-4), label = c("local mean"), color=("black"), size=notationsize)

#merge plots for global and local mean spectra
main_plot <- cowplot::plot_grid(
  full_plot2 + theme(legend.position="none", axis.title.x = element_blank(),  axis.text.x = element_blank()),
  full_plot1 + theme(legend.position="none"),
  nrow = 2, 
  rel_heights = c(0.45,0.55)
)

#add legend to merged plots
plot_main <- cowplot::plot_grid(
  main_plot,
   plot_grid(
      get_legend(for_model_legend1), 
      get_legend(for_obs_legend1),
      nrow = 2,
      align = "v"
    ), 
  nrow = 1, 
  rel_widths = c(0.8, 0.2)
)

#show plot
print(plot_main)
