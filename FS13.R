source("helpers/init.R")
source("helpers/functions.R")
meta.res_tbb <-readRDS("helpers/meta_res.Rds")
GMST <- readRDS("data/GMST_tbb.Rds")

split <- function(tibble){
    tibble <-  rowid_to_column(tibble, "ID")
    tbb1 <- tibble %>% filter(signal %in% c(signal_tbb %>% filter(type=="obs"&signal!="pages2k") %>% select(signal))$signal|grepl("highres",signal))
    tbb2 <- tibble %>% filter(!signal %in% tbb1$signal) %>% add_column(split="early") %>%
      mutate(data = purrr::map2(data, cut_time, function(.data, .cut_time){filter(.data, time <= .cut_time)})) 
    tbb3 <- tibble %>% filter(!signal %in% tbb1$signal) %>% add_column(split="late") %>%
      mutate(data = purrr::map2(data, cut_time, function(.data, .cut_time){filter(.data, time > .cut_time)})) 
    res <- rbind(tbb2, tbb3) %>% arrange(ID) %>% select(-ID) 
    return(res)
}

GMST <- inner_join(GMST, GMST %>% unnest(data) %>% group_by(signal) %>% summarise(cut_time=(min(time)+max(time))/2))

GMST_equi <- GMST %>% left_join(., meta.res_tbb) %>% mutate(interp.res = replace_na(interp.res, 1)) %>% 
  split %>% equidistant %>% select(signal, interp.res, EquiTS, cut_time, split)

GMST_spec <- tibble_spec(GMST_equi, 3, 2) %>% select(signal, interp.res, Spec, cut_time, split)

res <- tibble()
for(i in c("early", "late")){
  speclist <- tibble_to_list(GMST_spec %>% filter(split==i))
  names(speclist) <- unique(GMST_spec$signal)
  pages2k_spectra <- speclist[which(names(speclist) == "BHM"):which(names(speclist) == "PCR")]
  speclist[which(names(speclist) == "BHM"):which(names(speclist) == "PCR")] <- NULL
  speclist$pages2k <- MeanSpec(pages2k_spectra)$spec
  speclist <- lapply(speclist, function(x) AddConfInterval(x))
  speclist_smoothed <- lapply(speclist, function(x) LogSmooth(x, df=0.005, removeFirst=1, removeLast=floor(3*length(x$freq)/1000)))
  speclist_smoothed_tbb <- list_to_tibble(speclist_smoothed) %>% add_column(split=i)
  res <- rbind(res, speclist_smoothed_tbb)
}

res <- res %>% filter(!model %in% c("MPI-M_cont", "CESM_LM_cont")) %>% rename(., signal=model) %>%
  inner_join(., signal_tbb) %>% 
  arrange(., signal) %>% arrange(., type)

levels <- res %>% group_by(signal, type)  %>% select(-data) %>% filter(split=="early") %>% arrange(., order) 
cntlevels <- levels %>% group_by(type) %>% filter(split=="early") %>% count() 
cut_colors <-setNames(c(rev(brewer.pal(cntlevels[which(cntlevels$type=="model"),]$n+1, "YlGnBu"))[1:cntlevels[which(cntlevels$type=="model"),]$n], rev(brewer.pal(cntlevels[which(cntlevels$type=="obs"),]$n +1, "OrRd")[2:(cntlevels[which(cntlevels$type=="obs"),]$n+1)]), "grey70"), c(levels$signal, "mean"))

yrs.period <- rev(c(0.0001, 0.001, 0.01,  0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
yrs.labels <- rev(c(TeX('$10^{-4}$'),TeX('$10^{-3}$'), TeX('$10^{-2}$'), TeX('$10^{-1}$'), TeX('$10^{0}$'), TeX('$10^{1}$'), TeX('$10^{2}$'), TeX('$10^{3}$'), TeX('$10^{4}$'), TeX('$10^{5}$'), TeX('$10^{6}$')))
res %>% filter(!signal %in% c("HadCM3", "Trace21k_orb")) %>% unnest("data") %>%
        ggplot(aes(x = 1/freq, y = spec)) + 
        facet_wrap(~signal, nrow=2) +
        geom_ribbon(alpha=0.2, aes(ymin=lim.1, ymax=lim.2, fill=split), linetype = 0) + guides(fill=FALSE) +
        geom_line(aes(color=split), size=0.4) + theme_td() +
        guides(fill=FALSE)  +
        scale_fill_manual(values=c("#0F2080", "#F5793A")) + 
        scale_color_manual(values=c("#0F2080", "#F5793A"), name="selection:")  +
        scale_y_log10(name=TeX('PSD $S(\\tau)\\,(K^2 yr)$ '), label = scales::trans_format("log10", scales::math_format(10^.x)), expand=c(0.0, 0.0), limits=c(2e-4, 2e1), sec.axis = dup_axis(name = NULL, labels = NULL))  +
        scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, name=TeX('period $\\tau (yr)$'), expand=c(0.0, 0.0), limits=c(1e3, 1e-3), sec.axis = dup_axis(name = NULL, labels = NULL)) +
        theme(legend.key.size = unit(0.5, "cm"),
        legend.position = c(0.85, 0.25), 
          legend.title = element_text("selection"))
