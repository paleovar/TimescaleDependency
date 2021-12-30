source("helpers/init.R")
source("helpers/functions.R")

#--------------------------------------------------------------#
#FS6
#load and process data for plotting data
##local mean spectra
rmst_spec <- readRDS("data/local_mean_spectra.Rds")
rmst_spec <- rmst_spec %>% unnest(data) %>% filter(between(freq, 1/500, 1/1.5)) %>% group_by(signal) %>% nest() %>% 
inner_join(.,signal_tbb) %>% filter(signal %in% c("CESM_LM", "CESM_LM_cont", "MPI-M","MPI-M_cont"))  %>% ungroup() %>% select(-signal) %>% rename(signal=alt_name)
##global mean spectra
tmp_spec <- readRDS("data/global_mean_spectra.Rds")
tmp_spec <- tmp_spec %>% unnest(data) %>% filter(between(freq, 1/500, 1/1.5)) %>% group_by(signal) %>% nest() %>% inner_join(.,signal_tbb)

#define plotting parameters
levels <- tmp_spec %>% arrange(., order) %>% filter(!signal %in% c("CESM_LM_cont", "MPI-M_cont")) 
cntlevels <- levels %>% group_by(type) %>% count()
cut_colors <-setNames(c(rev(brewer.pal(cntlevels[which(cntlevels$type=="model"),]$n, "YlGnBu"))[1:cntlevels[which(cntlevels$type=="model"),]$n], rev(brewer.pal(cntlevels[which(cntlevels$type=="obs"),]$n +2, "OrRd")[2:(cntlevels[which(cntlevels$type=="obs"),]$n+2)])), c(levels$alt_name))
cut_colors <- setNames(c(rep(cut_colors[["CESM-LME 1"]],2), rep(cut_colors[["MPI-M LM"]], 2)), 
                       c("CESM-LME 1", "CESM-LME 1 (cntl)", "MPI-M LM", "MPI-M LM (cntl)"))
cut_colors <- setNames(c(rep("#0F2080",2), rep("#85C0F9", 2)), 
                       c("CESM-LME 1", "CESM-LME 1 (cntl)", "MPI-M LM", "MPI-M LM (cntl)"))
#create plot
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

#--------------------------------------------------------------#
#FS7
#load and process data for plotting data
tmp_spec <- readRDS("data/supp/CESM_spectra.Rds") %>% rename(signal=model)

#create plot
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

#--------------------------------------------------------------#
#FS8
#load required packages and data, and define functions
library(purrr)
prxlist <- readRDS("data/proxylist.Rds")
prxspec <- readRDS("data/proxy_spectra.Rds") 
prxscaling <- readRDS("data/beta_N100.Rds") %>% filter(signal=="pages2k")
tscales <- c(10,200)
stot <- function(specmtm){
  if("dof" %in% names(specmtm)){
    specmtm <- AddConfInterval(specmtm)
  return(tibble(spec=specmtm$spec, 
         freq=specmtm$freq,
         dof=specmtm$dof,
         lim.1=specmtm$lim.1, 
         lim.2=specmtm$lim.2))
    }
  if(!"dof" %in% names(specmtm)){
    return(tibble(spec=specmtm$spec, 
                  freq=specmtm$freq))
    }
}

#get timesteps of temporal resolution from irregulary sample proxy data
diffs.norm <- prxlist %>% mutate(diffs=purrr::map(data, function(x) diff(index(x)))) %>% select(-data) %>% rename(data=diffs)
tmp <- inner_join(prxspec, inner_join(diffs.norm, prxscaling, by="Name"), by="Name")

#create surrogate timeseries with powerlaw scaling and temporal resolution as proxy records
tbb <- list()
tbb <- tmp %>% mutate(tmp = purrr:::map2(data, slope, function(.data, .slope) sample_from_proxy(.data, (-1)*.slope))) %>%
  mutate(powerlaw = purrr::map(tmp, function(.tmp) .tmp[["powerlaw"]]), 
         sample = purrr::map(tmp, function(.tmp) .tmp[["sample"]]))   
#interpolate the time series
tbb <- tbb %>% mutate(ts = map2(sample, interp.res, function(.sample, .interp.res) 
  MakeEquidistant(index(.sample), coredata(.sample), dt=.interp.res))) 
tbb <- tbb %>%  mutate(comparison = map2(powerlaw, interp.res, function(.powerlaw, .interp.res)
    MakeEquidistant(index(as.ts(.powerlaw)), coredata(as.ts(.powerlaw)), dt = .interp.res))) 
#compute spectrum for block-averaged, equidistant artificial proxy series (ts) and for the comparison   
tbb <- tbb %>% mutate(spec = purrr::map(ts, function(.ts) SpecMTM(na.approx(.ts), k=3, nw=2, detrend=TRUE))) 
tbb <- tbb %>% mutate(comparison.spec =  map2(comparison, ts, function(.comparison, .none){ SpecMTM(.comparison, k=3, nw=2, detrend=TRUE)}))

#compute scaling coefficients
tbb <- tbb %>% mutate(temp = map2(spec, comparison.spec, function(.spec, .none) 
  SlopeFit(.spec, 1/200, 1/10))) %>%    
  mutate(simslope = map2(temp, comparison.spec, function(.temp, .none) .temp$slope)) %>%     
  mutate(slopesd = map2(temp, comparison.spec, function(.temp, .none) .temp$slopesd)) %>% 
  mutate(compslope = map2(comparison.spec, temp, function(.comparison.spec, .none)
    SlopeFit(.comparison.spec, 1/200, 1/10)$slope)) %>% 
  select(-temp, -tmp)

#create example plot
##parameters
ex <- which(tbb$Name == "NAm-LakeMina")
yrs.period <- rev(c(0.0001, 0.001, 0.01,  0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
yrs.labels <- rev(c(TeX('$10^{-4}$'),TeX('$10^{-3}$'), TeX('$10^{-2}$'), TeX('$10^{-1}$'), TeX('$10^{0}$'), TeX('$10^{1}$'), TeX('$10^{2}$'), TeX('$10^{3}$'), TeX('$10^{4}$'), TeX('$10^{5}$'), TeX('$10^{6}$')))  
colors <- c("surrogate" = "#0F2080", "block averaged surrogate" = "#85C0F9", "raw signal" = "#F5793A")
colors <- c(colors, "smoothed spectrum"="#0F2080", "fit"="black")
linetype <- c("surrogate" = "dotted", "block averaged surrogate" = "solid", "raw signal" = "dashed")
linetype <- c(linetype, "smoothed spectrum"="dotted", "fit"="longdash", "raw signal" = "dashed")
pointshape <- c("surrogate" = 17,"block averaged surrogate" = 16, "raw signal" = 15)
##processing before plotting
raw <- SpecMTM(na.remove(MakeEquidistant(index(c(prxlist %>% filter(Name==tbb$Name[[ex]]))$data[[1]]), coredata(c(prxlist %>% filter(Name==tbb$Name[[ex]]))$data[[1]]), dt=tbb$interp.res[[ex]])))
bavg_surrogate <- tbb$comparison.spec[[ex]]
bavg_surrogate_fit <- SlopeFit(bavg_surrogate, 1/200, 1/10, bDebug=FALSE)
slope <-bavg_surrogate_fit$slope
slopesd <-bavg_surrogate_fit$slopesd
int <- bavg_surrogate_fit$intercept
temp <- paste(paste('beta[fit] == ', round(-1*slope,2)), " %+-% ", round(slopesd,2))

#create plot
ggplot() + theme_td() +
  geom_ribbon(data=stot(raw), alpha=0.1, aes(x = 1/freq, ymin=lim.1, ymax=lim.2,  fill="raw signal")) +
  geom_ribbon(data=stot(bavg_surrogate), alpha=0.2, aes(x = 1/freq, ymin=lim.1, ymax=lim.2,  fill="block averaged surrogate")) +
  scale_fill_manual(values=colors) +
  scale_y_log10(TeX('PSD $S(\\tau)\\, (K^2 yr)$ '), label = trans_format("log10", math_format(10^.x)),expand=c(0., 0.1), sec.axis = dup_axis(name = NULL, labels = NULL))  + 
  scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels,  expand=c(0., 0.), sec.axis = dup_axis(name = NULL, labels = NULL), name=TeX('period $\\tau\\,(yr)$')) + 
  geom_line(data=stot(raw), aes(x=1/freq,y=spec, color="raw signal", linetype="raw signal"), size=pointsize) +
  geom_line(data=stot(bavg_surrogate), aes(x=1/freq, y=spec, color='block averaged surrogate', linetype="block averaged surrogate"), alpha=0.8, size=pointsize) +
  geom_line(data=stot(bavg_surrogate_fit), aes(x=1/freq, y=spec, color='smoothed spectrum', linetype="smoothed spectrum"), alpha=0.8, size=pointsize) +
  geom_line(data= stot(bavg_surrogate), aes(x=1/freq, y=((1/freq)**(-slope))*(exp(int)), linetype="fit"), color="black") +
  geom_vline(xintercept = c(200, 10), color="grey50", alpha=0.5) + 
  theme(legend.position=c(0.41, 0.34)) +
  scale_color_manual(values = colors, breaks=c("raw signal", "block averaged surrogate", "smoothed spectrum"), labels=c("raw spectrum", "surrogate spectrum", "smooothed surrogate spectrum"))  +
  scale_linetype_manual(values=linetype) +
  annotate(geom="text", x=50, y=0.003, label=expression(paste("Proxy: NAm-LakeMina, ", beta, "=1.31,", Delta, "t=3.98")), size=notationsize) +
  annotate(geom="text", x=28, y=100, label= temp,  size=notationsize, parse=TRUE) +
  guides(colour = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dotted")))) +
  guides(fill=FALSE, linetype=FALSE)

#--------------------------------------------------------------#
#FS9
#load data
tmp_spec <- readRDS("data/supp/pages_spectra_selection.Rds") %>% rename(signal=model)

#create plot
tmp_spec %>% plot_spec(ylims=c(0.05, 30), xlims=c(1.e3, 2), name.col="cutoff", name.fill="res", name.line="res", name.alpha="cutoff") +
    scale_fill_manual(values=c("#F5793A", "#0F2080")) +
    scale_linetype_manual(values=c("dashed", "solid"), name="coverage:", labels=c("PI", "hist")) +
    scale_color_manual(values=c( "#F5793A",  "#0F2080"), name="selection:") +   
    scale_alpha_manual(values=rep(0.1, 4)) +
    guides(alpha=F, linetype = guide_legend(order = 2), color = guide_legend(order = 1)) +
    theme(legend.key.size = unit(0.5, "cm"),
          legend.title = element_text("selection"),
          legend.position = c(0.3, 0.15))
rm(tmp_spec)
gc()

#--------------------------------------------------------------#
#FS10
#load required packages and data
library(ggridges)
library(forcats)
tbb <- rbind(
    readRDS("data/beta_N100.Rds"), 
    readRDS("data/supp/beta_N100_1850.Rds") 
)

#get two data sets with and without recent global warming
one <- tbb %>% filter(cutoff ==2020) %>% rename(slope2020=slope, slopesd_new2020=slopesd_new) %>% select(-slopesd)  %>% select(-cutoff)
two <- tbb %>% filter(cutoff ==1850) %>% select(-cutoff, ID)

#get difference in \beta
diffs_warming <- inner_join(one, two, by=c("Name", "Archive", "long", "lat", "signal", "ID", "scale"))  %>% mutate(slope_diff=-1*(slope2020-slope)) %>% filter(!signal%in%c("MPI-M_cont", "CESM_LM_cont")) 

#prepare plot
levels <- diffs_warming %>% group_by(signal) %>% nest() %>% select(signal) %>% inner_join(.,signal_tbb) %>% arrange(., order) 
cntlevels <- levels %>% group_by(type) %>% count()
cut_colors <-setNames(c(rev(brewer.pal(cntlevels[which(cntlevels$type=="model"),]$n, "YlGnBu"))[1:cntlevels[which(cntlevels$type=="model"),]$n], rev(brewer.pal(cntlevels[which(cntlevels$type=="obs"),]$n +2, "OrRd")[2:(cntlevels[which(cntlevels$type=="obs"),]$n+2)])), c(levels$alt_name))
diffs_warming <- diffs_warming %>% filter(!signal%in%c("Trace21k_orb", "HadCM3")) %>% inner_join(signal_tbb, by="signal") %>% mutate(signal = factor(signal, levels = unique(.$signal))) %>% mutate(signal = fct_rev(as_factor(signal)))

#create plot
ggplot(diffs_warming, aes(x=slope_diff, fill= fct_reorder(diffs_warming$alt_name, diffs_warming$order, max), y= fct_reorder(diffs_warming$alt_name, diffs_warming$order, max))) + 
  geom_density_ridges(scale = 1.) + 
  xlab(TeX("$\\beta_{\\mathrm{hist}} - \\beta_{\\mathrm{PI}}$")) +
  scale_y_discrete(limits=rev(c("pages2k", unique((levels %>% filter(!signal %in% c("pages2k", "Trace21k_orb", "HadCM3")))$alt_name))), 
                              breaks=c("pages2k", unique((levels %>% filter(!signal %in% c("pages2k", "Trace21k_orb", "HadCM3")))$alt_name)), 
                   labels=c("Proxy", unique((levels %>% filter(!signal%in% c("pages2k", "Trace21k_orb", "HadCM3")))$alt_name)), expand=c(0.12,0.1)) +
  theme_td() +
  theme(panel.background = element_rect(fill = "transparent", colour = "black", size=.8)) +
  scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(legend.position="none",
        axis.title.y = element_blank()) +
  scale_fill_manual(values=cut_colors) 

