source("helpers/init.R")
library(cowplot)
library(forcats)

gmspectra <- readRDS("data/global_mean_spectra.Rds")
lmspectra <- readRDS("data/local_mean_spectra.Rds")

mean_df <- gmspectra %>% filter(signal == "mean") %>% unnest(data)
gmspectra <- gmspectra %>% filter(!signal %in% c("MPI-M_cont", "CESM_LM_cont"))  %>%
  inner_join(., signal_tbb) %>% 
  arrange(., signal) %>% arrange(., type)

levels <- gmspectra %>% group_by(signal, type)  %>% select(-data)  %>% arrange(., order)
cntlevels <- levels %>% group_by(type) %>% count()
cut_colors <-setNames(c(rev(brewer.pal(cntlevels[which(cntlevels$type=="model"),]$n+1, "YlGnBu"))[1:cntlevels[which(cntlevels$type=="model"),]$n], rev(brewer.pal(cntlevels[which(cntlevels$type=="obs"),]$n +1, "OrRd")[2:(cntlevels[which(cntlevels$type=="obs"),]$n+1)]), "black"), c(levels$signal, "mean"))

full_plot1 <- gmspectra %>% unnest(data) %>% 
  ggplot() + theme_td() +
  geom_abline(intercept = c(-5, -0.5), slope = c(-1, -2), lty=2, size=0.2) +
  geom_ribbon(alpha=0.1, aes(x = 1/freq, ymin=lim.1, ymax=lim.2, fill=signal)) +
  geom_ribbon(data=mean_df, alpha=0.1, aes(x = 1/freq, ymin=lim.1, ymax=lim.2, fill="mean")) +
  geom_line(size=pointsize, aes(x = 1/freq, y = spec, color=signal))  +
  scale_color_manual(values = cut_colors) +
  guides(fill=FALSE) +
  geom_line(data=mean_df, aes(x=1/freq,y=spec), color='black',size=pointsize) +
  scale_fill_manual(values=cut_colors) + 
  theme(legend.position="none") +
  annotate("text", x = c(0.02, 20), y=c(0.01, 0.00003), label = c("β=2", "β=1"), size=notationsize) +
  annotate("text", x = c(350), y=c(1e-7), label = c("global mean"), color=("black"), size=notationsize) +
  annotate("text",x = c(500), y=c(500), label = c("(b)"), size=notationsize)  + theme_td() +
  scale_y_log10(name=TeX('PSD S($\\tau$) ($K^2 yr$)'), breaks = c(0.000001, 0.0001, 0.01, 1, 100), labels = c(TeX('$10^{-6}$'), TeX('$10^{-4}$'), TeX('$10^{-2}$'), TeX('$10^{0}$'), TeX('$10^{2}$')), expand=c(0.0, 0.0), limits=c(1e-8, 2e3), sec.axis = dup_axis(name = NULL, labels = NULL))  +
  scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, name=TeX("period $\\tau$ ($yr$)"), expand=c(0.0, 0.), limits=c(1e3, 1e-3),sec.axis = dup_axis(name = NULL, labels = NULL)) 


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

for_obs_legend1 <- gmspectra %>% unnest(data) %>%
  filter(signal %in% obs) %>%
  ggplot(aes(x = 1/freq, y = spec, ymin=lim.1, ymax=lim.2,  fill=signal)) + 
  geom_line(aes(color=signal)) +
  scale_color_manual(values = cut_colors, name = "observation-based:") + 
  theme(legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.65, 'lines'),
        legend.title = element_text(size=9*0.75, face="bold"),
        legend.text = element_text(size=9*0.75))

lmspectra_mean_df <- lmspectra %>% filter(signal == "mean") %>% unnest(data)
lmspectra <- lmspectra %>% filter(!signal %in% c("MPI-M_cont", "CESM_LM_cont")) %>% inner_join(., signal_tbb) %>% arrange(., signal) %>% arrange(., type)

levels <- lmspectra %>% group_by(signal, type) %>% select(-data) %>% arrange(., order)
cntlevels <- levels %>% group_by(type) %>% count()
cut_colors <-setNames(c(rev(brewer.pal(cntlevels[which(cntlevels$type=="model"),]$n+1, "YlGnBu"))[1:cntlevels[which(cntlevels$type=="model"),]$n], rev(brewer.pal(cntlevels[which(cntlevels$type=="obs"),]$n +1, "OrRd")[2:(cntlevels[which(cntlevels$type=="obs"),]$n+1)]), "black"), c(levels$signal, "mean"))

full_plot2 <- lmspectra %>% unnest(data) %>% 
  ggplot() + theme_td() +
  geom_abline(intercept = c(-3, 2.5), slope = c(-1, -2), lty=2, size=0.2) +
  scale_y_log10(name=TeX('PSD S($\\tau$) ($K^2 yr$)'),breaks = c(0.000001, 0.0001, 0.01, 1, 100), labels = c(TeX('$10^{-6}$'), TeX('$10^{-4}$'), TeX('$10^{-2}$'), TeX('$10^{0}$'), TeX('$10^{2}$')), limits=c(1e-5, 2e3), expand=c(0.0, 0.0), sec.axis = dup_axis(name = NULL, labels = NULL))  +
  scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, limits=c(1e3, 1e-3), expand=c(0.0, 0.), name=TeX('\\textbf{period} / $yrs$'), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  geom_ribbon(alpha=0.3, aes(x = 1/freq, ymin=lim.1, ymax=lim.2, fill=signal)) +
  geom_ribbon(data=lmspectra_mean_df, alpha=0.3, aes(x = 1/freq, ymin=lim.1, ymax=lim.2, fill="mean")) +
  geom_line(size=pointsize, aes(x = 1/freq, y = spec, color=signal))  +
  scale_color_manual(values = cut_colors) +
  guides(fill=FALSE) +
  geom_line(data=lmspectra_mean_df, aes(x=1/freq,y=spec), color='black',size=pointsize) +
  scale_fill_manual(values=cut_colors) + #
  annotate("text", x = c(0.02, 5), y=c(10, 0.0003), label = c("β=2", "β=1"), size=notationsize) +
  annotate("text",x = c(500), y=c(500), label = c("(a)"), size=notationsize) +
  annotate("text",x = c(350), y=c(1e-4), label = c("local mean"), color=("black"), size=notationsize) +
  theme(legend.position="none") + theme_td()

main_plot <- cowplot::plot_grid(
  full_plot2 + theme(legend.position="none", axis.title.x = element_blank(),  axis.text.x = element_blank()),
  full_plot1 + theme(legend.position="none"),
  nrow = 2, 
  rel_heights = c(0.45,0.55)
)

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

print(plot_main)
