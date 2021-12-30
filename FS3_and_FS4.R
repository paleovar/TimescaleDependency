source("helpers/init.R")
source("helpers/functions.R")
#--------------------------------------------------------------#
#FS3

tbb <- readRDS("data/pages_meta_scaling.Rds")
pages.prxlist <- readRDS("data/proxylist.Rds") %>% inner_join(., tbb,  by="Name") 
#process data for plotting
pages.prxlist$data2 <- lapply(pages.prxlist$data, function(x) tibble(time=index(x) / 1000, val=coredata(x) - mean(na.approx(coredata(x)))))
dat <- pages.prxlist %>% unnest(data2)
rm(pages.prxlist)

#plot temperature anomalys
ggplot(dat, aes(x=time, y=val)) + 
    geom_line(size=0.1) + xlab("time (kyr CE)") + ylab(TeX("Temperature anomaly (K)")) + geom_line(size=0.1) +
    scale_x_continuous( labels = scales::number_format(accuracy = 0.1), expand=c(0.1,0.1)) +
    scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
    theme_td() + 
    theme(legend.position = "none") + facet_wrap(~ID, scales = "free_x") 

#--------------------------------------------------------------#
#FS4
#load data
tmp_spec <- readRDS("data/proxy_spectra.Rds") 

#plotting parameters
yrs.period <- rev(c(0.0001, 0.001, 0.01,  0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
yrs.labels <- rev(c(TeX('$10^{-4}$'),TeX('$10^{-3}$'), TeX('$10^{-2}$'), TeX('$10^{-1}$'), TeX('$10^{0}$'), TeX('$10^{1}$'), TeX('$10^{2}$'), TeX('$10^{3}$'), TeX('$10^{4}$'), TeX('$10^{5}$'), TeX('$10^{6}$')))

#plot spectra
ggplot(tmp_spec %>% unnest(Spec)) + 
    facet_wrap(~ID) +
    geom_line(aes(x=1/freq, y=temp),color="black", size=0.1) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)), name=TeX('PSD $S(\\tau)\\, (K^2 yr)$ '), expand=c(0.2, 0.2), 
                  sec.axis = dup_axis(name = NULL, labels = NULL))  + 
    scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels,   
                       expand=c(0.05, 0.05), limits=c(1e3, 3), name=TeX('period $\\tau\\,(yr)$')) + 
    theme(legend.position="none")  +
    annotate(geom="rect", xmin = c(10), ymin = c(0), xmax = c(200), ymax = c(Inf), fill="#f5793a", alpha=0.3)+ 
    theme_td()
