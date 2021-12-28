#FS6

source("helpers/init.R")
source("helpers/functions.R")
rmst_spec <- readRDS("data/local_mean_spectra.Rds")
rmst_spec <- rmst_spec %>% unnest(data) %>% filter(between(freq, 1/500, 1/1.5)) %>% group_by(signal) %>% nest() %>% 
inner_join(.,signal_tbb) %>% filter(signal %in% c("CESM_LM", "CESM_LM_cont", "MPI-M","MPI-M_cont"))  %>% ungroup() %>% select(-signal) %>% rename(signal=alt_name)

tmp_spec <- readRDS("data/global_mean_spectra.Rds")
tmp_spec <- tmp_spec %>% unnest(data) %>% filter(between(freq, 1/500, 1/1.5)) %>% group_by(signal) %>% nest() %>% inner_join(.,signal_tbb)

levels <- tmp_spec %>% arrange(., order) %>% filter(!signal %in% c("CESM_LM_cont", "MPI-M_cont")) 

cntlevels <- levels %>% group_by(type) %>% count()
cut_colors <-setNames(c(rev(brewer.pal(cntlevels[which(cntlevels$type=="model"),]$n, "YlGnBu"))[1:cntlevels[which(cntlevels$type=="model"),]$n], rev(brewer.pal(cntlevels[which(cntlevels$type=="obs"),]$n +2, "OrRd")[2:(cntlevels[which(cntlevels$type=="obs"),]$n+2)])), c(levels$alt_name))
cut_colors <- setNames(c(rep(cut_colors[["CESM-LME 1"]],2), rep(cut_colors[["MPI-M LM"]], 2)), 
                       c("CESM-LME 1", "CESM-LME 1 (cntl)", "MPI-M LM", "MPI-M LM (cntl)"))
cut_colors <- setNames(c(rep("#0F2080",2), rep("#85C0F9", 2)), 
                       c("CESM-LME 1", "CESM-LME 1 (cntl)", "MPI-M LM", "MPI-M LM (cntl)"))

tmp_spec %>% ungroup() %>% select(-signal) %>% rename(signal=alt_name) %>% filter(signal%in% c("CESM-LME 1", "CESM-LME 1 (cntl)", "MPI-M LM", "MPI-M LM (cntl)")) %>%
    plot_spec(ylims=c(3e-4,8), xlims=c(1000, 0.5), name.line="signal") + guides(fill=FALSE, alpha=FALSE, linetype=FALSE) +
    theme(legend.position = c(0.3, 0.2)) +
    geom_ribbon(data=rmst_spec %>% unnest(data), aes(x = 1/freq, ymin=lim.1, ymax=lim.2, fill=signal), alpha=0.5) +
    geom_line(data=rmst_spec %>% unnest(data), aes(x = 1/freq, y = spec, color=signal, linetype=signal), size=pointsize) +
    scale_fill_manual(values=cut_colors) +
    scale_color_manual(values=cut_colors) +
    scale_linetype_manual(values=c("CESM-LME 1"="solid", "CESM-LME 1 (cntl)"="dashed", "MPI-M LM"="solid", "MPI-M LM (cntl)"="dashed")) +
    annotate(geom="text", x=c(8,8), y=c(0.006, 2), label=c("global", "local"), size=9/2.8, col="black") +
    guides(colour = guide_legend(override.aes = list(linetype = c("solid", "dashed", "solid", "dashed"))))
rm(tmp_spec)
gc()
#---------------------#
#FS7
tmp_spec <- readRDS("data/CESM_spectra.Rds") %>% rename(signal=model)
tmp_spec %>% plot_spec(ylims=c(0.00002, 2000), xlims=rev(c(4e2,1.2e-3)), name.line="signal") +
guides(fill=FALSE, alpha=FALSE, linetype=FALSE) +
scale_linetype_manual(values=c("CESM_LM"="solid",  "CESM_LM_highres"="solid", "CESM_mean"="dashed")) +
theme(legend.position = c(0.3, 0.15)) +
scale_color_manual(values=c("CESM_LM"="#0F2080",  "CESM_LM_highres"="#85C0F9", "CESM_mean"="#F5793A"), 
                     labels=c("CESM-LME 1",  "CESM-LME 1 (hr)", "mean spectrum"),
                     breaks=c("CESM_LM",   "CESM_LM_highres", "CESM_mean")) +
guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed"))))
rm(tmp_spec)
gc()
#---------------------#
#FS9
tmp_spec <- readRDS("data/pages_spectra_selection.Rds") %>% rename(signal=model)

tmp_spec %>% plot_spec(ylims=c(0.05, 30), xlims=c(1.e3, 2), name.col="cutoff", name.fill="res", name.line="res", name.alpha="cutoff") +
    scale_fill_manual(values=c("#F5793A", "#0F2080")) +
    scale_linetype_manual(values=c("dashed", "solid"), name="coverage:", labels=c("PI", "hist")) +
    scale_color_manual(values=c( "#F5793A",  "#0F2080"), name="selection:") +   
    scale_alpha_manual(values=rep(0.1, 4)) +
    guides(alpha=F, linetype = guide_legend(order = 2), color = guide_legend(order = 1)) +
    theme(legend.key.size = unit(0.5, "cm"),
          legend.title = element_text("selection"),
          legend.position = c(0.3, 0.15))
