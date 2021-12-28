source("helpers/init.R")
source("helpers/functions.R")
tbb <- readRDS("data/pages_meta_scaling.Rds")

#pages.prxlist <- readRDS(paste0(pages.dir, "/pages.prxlist_", min.res, "_", min.range, "_", max.hiat, ".Rds"))
pages.prxlist <- readRDS("data/proxylist.Rds") %>% inner_join(., tbb,  by="Name") 

pages.prxlist$data2 <- lapply(pages.prxlist$data, function(x) tibble(time=index(x) / 1000, val=coredata(x) - mean(na.approx(coredata(x)))))
  
dat <- pages.prxlist %>% unnest(data2)

ggplot(dat, aes(x=time, y=val)) + 
    geom_line(size=0.1) + xlab("time (kyr CE)") + ylab(TeX("Temperature anomaly (K)")) + geom_line(size=0.1) +
    scale_x_continuous( labels = scales::number_format(accuracy = 0.1), expand=c(0.1,0.1)) +
    scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
    theme_td() + 
    theme(legend.position = "none") + facet_wrap(~ID, scales = "free_x") 
