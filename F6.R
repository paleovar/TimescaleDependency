source("helpers/init.R")
source("helpers/functions.R")

l <- readRDS("data/transfer.Rds")

plot_spec(l, ylims=c(0.001, 1), xlims=c(370, 2.1))

plot_spec <- function(tibble, ylims, xlims, name.y=TeX('PSD $S(\\tau)\\,(K^2 yr)$ '), name.x=TeX('period $\\tau (yr)$'), name.data = "data"){
    yrs.period <- rev(c(0.0001, 0.001, 0.01,  0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
    yrs.labels <- rev(c(TeX('$10^{-4}$'),TeX('$10^{-3}$'), TeX('$10^{-2}$'), TeX('$10^{-1}$'), TeX('$10^{0}$'), TeX('$10^{1}$'), TeX('$10^{2}$'), TeX('$10^{3}$'), TeX('$10^{4}$'), TeX('$10^{5}$'), TeX('$10^{6}$')))
    tibble %>% unnest(name.data) %>%
        ggplot(aes(x = 1/freq, y = spec, color=signal)) +
        geom_ribbon(alpha=0.2, aes(ymin=lim.1, ymax=lim.2, fill=signal), linetype = 0) + guides(fill=FALSE) +
        geom_line(size=0.4) + theme_td() +
        scale_y_log10(name=name.y, label = scales::trans_format("log10", scales::math_format(10^.x)), expand=c(0.0, 0.0), limits=ylims, sec.axis = dup_axis(name = NULL, labels = NULL))  +
        scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, name=name.x,expand=c(0.0, 0.0), limits=xlims, sec.axis = dup_axis(name = NULL, labels = NULL)) 
    }

sample_gain <- l %>% rename(lim.1=sd_down, lim.2=sd_up, signal=name, spec=m) %>% group_by(signal) %>% nest() %>% 
  plot_spec(ylims=c(0.001, 1), xlims=c(370, 2.1),name.y=TeX('gain G^2 ($\\tau$) ($K^2/(W^2 m^{-4})$)')) + 
          geom_hline(aes(x=1/freq, yintercept=var, color=signal), linetype="dashed") +
          scale_color_manual(values=c("grey70",  "#6699cc", "#ffcc65"), labels=c("models w/ ENSO",  "models w/o ENSO", "reconstructions"), 
                     breaks=c("models_wENSO", "models_woENSO", "recons")) + 
          scale_fill_manual(values=c("grey50",  "#6699cc", "#ffcc65")) + 
          annotate("text",x = c(280), y=c(0.9), label = c("(b)"), size=notationsize) +
          theme(legend.position=c(0.27, 0.18)) + 
          guides(color=guide_legend(override.aes=list(fill=NA)))

print(sample_gain)

l_raw <- readRDS("data/sample_spec.Rds")

sample_raw <- l_raw %>% rename(lim.1=sd_down, lim.2=sd_up, signal=name, spec=m) %>% group_by(signal) %>% nest() %>%
 plot_spec(xlims=c(370, 2.1), ylims=c(5e-4, 50)) + 
  scale_color_manual(values=c("#009966", "grey70",  "#6699cc", "#ffcc65"), labels=c(unname(TeX("joint forcing")), unname(TeX("model simulations $(M_0, M_+)$")), unname(TeX("model simulations $(M_0)$")),  unname(TeX("observation-based"))), 
                     breaks=c("forc", "models_wENSO", "models_woENSO", "recons")) + # "#6699cc", "
  scale_fill_manual(values=c( "#009560", "grey50",  "#6699cc", "#ffcc65")) +
  scale_y_log10(name=TeX('PSD $S_T (\\tau)$ ($K^2 yr$)'), breaks=rev(c(0.0001, 0.001, 0.01,  0.1, 1, 10, 100, 1000, 10000, 100000, 1000000)),
    labels=rev(c(TeX('$10^{-4}$'),TeX('$10^{-3}$'), TeX('$10^{-2}$'), TeX('$10^{-1}$'), TeX('$10^{0}$'), TeX('$10^{1}$'), TeX('$10^{2}$'), TeX('$10^{3}$'), TeX('$10^{4}$'), TeX('$10^{5}$'), TeX('$10^{6}$'))),
    expand=c(0.01,0.4),  
    sec.axis = dup_axis(name = TeX('PSD $S_F (\\tau)$ (W^2 m^{-4} yr$)'), labels = NULL)) +
  theme(legend.text.align = 0,
    axis.title.y.right = element_text(color = "#009560"),
    legend.position=c(0.2,0.2) 
  ) +
  annotate("text",x = c(280), y=c(15), label = c("(a)"), size=notationsize) +
  annotate("text",x = c(10), y=c(8), label = expression(S[F]), size=notationsize) +
  annotate("text",x = c(10), y=c(0.1), label = expression(S[T]), size=notationsize) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

print(sample_raw)

cowplot::plot_grid(
  sample_raw + 
  theme(axis.text.x = element_blank(), axis.title.x=element_blank()),
  sample_gain + theme(legend.position="none"),
  nrow = 2, 
  align="v",
  rel_heights = c(0.47, 0.53)
)
