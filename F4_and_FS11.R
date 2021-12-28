library(irr)
source("helpers/init.R")

#F4
tbb <- readRDS("data/beta_N100.Rds") 

kappa.stat <- function(target1, target2, method){
  tmp <- categories %>% select(Name, signal, category) %>% filter(signal %in% c(target1, target2)) %>% group_by(signal) %>% nest()
  if(any(tmp$data[[1]]$Name != tmp$data[[2]]$Name)){stop("Names not matching for cases studied")}
  if(method=="agreement"){
    rating <- data.frame(
      rtr1 = tmp$data[[1]]$category,
      rtr2 =tmp$data[[2]]$category
    )
    return(agree(rating)$value/100)
  }
  if(method=="kripp"){
    rating <- cbind(tmp$data[[1]]$category, tmp$data[[2]]$category)
    return(kripp.alpha(rating)$value)
  }
  if(method=="kappa"){
    rating <- cbind(tmp$data[[1]]$category, tmp$data[[2]]$category)
    return(kappa2(rating)$value)
  }
}

overlap <- function(mean, dev, comp.idx=1){
  dev <- abs(dev)
  tar <- cbind(mean-dev, mean+dev)
  comp <- tar[comp.idx, ]
  if(length(dim(tar))==0){
    return((max(tar) >= max(comp) && min(tar) <= max(comp)) || max(tar) <= max(comp) && max(tar) >= min(comp))
  }
  if(length(dim(tar))!=0){
    return(
      apply(tar, 1, function(x) ((max(x) >= max(comp) && min(x) <= max(comp)) || max(x) <= max(comp) && max(x) >= min(comp)))
    )
  }
}

tbb <- tbb %>% inner_join(signal_tbb %>% select(signal, alt_name), by="signal") %>% 
  add_column(alt_name2= factor(.$alt_name, as.character(signal_tbb$alt_name))) 

assign <- tbb %>% select(ID, long, lat, signal) %>% filter(signal =="pages2k") %>% arrange(.,desc(lat))

tbb$lat <- factor(tbb$lat, levels = assign$lat)
 
t <- c()
for(i in 1:length(assign$lat)){
  t <- rbind(t, paste0(assign$ID[[i]], " (", round(assign$long[[i]],1), ":", round(assign$lat[[i]],1), ")"))
}

facet.labs <- t
names(facet.labs) <- levels(tbb$lat)

tbb <- tbb %>% inner_join(signal_tbb) %>% arrange(., order) 

ggplot(data = tbb,aes(x = alt_name2, y = -slope, ymin=-slope+slopesd_new, ymax=-slope-slopesd_new, shape=Archive)) +   
  annotate(geom="rect", xmin = c(-Inf), ymin = c(1.), xmax = c(Inf), ymax = c(Inf), fill="red3", alpha=0.4)+ 
  annotate(geom="rect", xmin =  c(-Inf), ymin = c(-0.5), xmax = c(Inf), ymax = c(0.5),  alpha=0.3, fill="white") + 
  geom_errorbar(width = 1, colour="black") +  
  geom_point(data=tbb, aes(size=Archive, fill=Archive))+ 
  scale_y_continuous(expand=c(0.05,0), limits=c(-1, 3)) +
  facet_wrap(. ~ lat, nrow=4, labeller=labeller(lat=facet.labs)) +  labs(x = NULL, y = TeX("scaling coefficient $\\beta$")) +   
  coord_flip(expand = TRUE) +  
  scale_x_discrete(limits=rev(c("pages2k", unique((tbb %>% filter(signal!="pages2k"))$alt_name))), breaks=c("pages2k", unique((tbb %>% filter(signal!="pages2k"))$alt_name)), labels=c("Proxy", unique((tbb %>% filter(signal!="pages2k"))$alt_name))) + 
  scale_size_manual(values=rep(1.2, 4), guide=FALSE) +   
  scale_shape_manual(values=c(24, 23, 22, 25), breaks=c("model","tree","lake/marine sediment","documents"), labels=c("model","tree","lake/marine\nsediment","documents")) + scale_fill_manual(values = rep("black",4)) +   
  theme_td() +  
  guides(fill=FALSE) +
  guides(shape = guide_legend(override.aes = list(fill = rep("black",4)))) +
  theme(legend.key.size = unit(0.75,"line"),        
        legend.position =c(0.92, 0.12),   
        axis.text.x = element_text(angle = 0, hjust=1, size=7),        
        legend.text=element_text(size=9),         
        legend.title=element_blank(), 
        legend.key.height=unit(0.7, "cm"),
        legend.key.width=unit(0.2, "cm")) 

#---------------------#
#FS11
  
ns <- unique(tbb$Name)
res <- tibble()
for(i in ns){
  loc_tbb <- tbb %>% filter(Name==i) 
  loc_tbb <- loc_tbb %>% add_column(agree=overlap(loc_tbb %>% select(slope), loc_tbb %>% select(slopesd_new), which(loc_tbb$signal=="pages2k")))
  res <- rbind(res, loc_tbb)
} 

percent_agree <- tibble(Name=character(), value=numeric(), method=character())
for(sig in unique(tbb$signal)){
  tmp <- tibble(Name=sig,
                value = (res %>% filter(signal==sig) %>% count(cnt= (agree==TRUE)) %>% filter(cnt=="TRUE") %>% select(n))$n / nrow(res %>% filter(signal==sig)), 
                method="percent agree wi margins")
  percent_agree <- rbind(percent_agree, tmp)
}

percent_agree$value <- round(percent_agree$value, 2)

mean_error <- as.numeric(tbb %>% filter(signal =="pages2k") %>% summarise(mean=mean(slopesd_new)) )

categories <- tbb %>% 
  mutate(category = case_when(slope <= -1-mean_error~ -2, 
                              slope >= -1-mean_error& slope <= -1+mean_error~ -1,
                              slope >= -1+mean_error ~ 0))
                              
res <- list()
tbb <- tibble(Name=character(), value=numeric(), method=character())
for(method in c("agreement", "kripp", "kappa")){ 
  #l <-list()
  for(n in unique(categories$signal)){
    if(n=="pages2k"){next}
    tmp <- tibble(Name=n,
           value = kappa.stat("pages2k", n, method), 
           method=method)
    tbb <- rbind(tbb, tmp)
  }
}

tbb$value <- round(tbb$value, 2)
tbb <- tbb %>% rbind(., percent_agree)

tbb <- tbb %>% rename(signal=Name) %>% inner_join(signal_tbb %>% select(signal, alt_name), by="signal") %>% 
  add_column(alt_name2= factor(.$alt_name, rev((signal_tbb %>% arrange(., order))$alt_name))) 

tbb <- tbb %>% filter(method!="kripp") %>% add_column(method2= factor(.$method, rev(c("percent agree wi margins", "agreement", "kappa"))))

ggplot(tbb %>% filter(!alt_name2=="pages2k", method!="kripp")) + 
  geom_bar(aes(y=value, x=alt_name2, fill=method2), width=1.2/1.5, stat="identity", position="dodge") + 
  scale_y_continuous(limits=c(-0.2, 1), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(),panel.background = element_rect(fill = "transparent", colour = "black", size=1.2)) +
   scale_fill_manual(values=c("percent agree wi margins"="#0f2080", "agreement"="#85c0f9", "kappa"="#f5793a"),
                    breaks=c("percent agree wi margins", "agreement", "kappa"), 
                    labels=c(bquote(paste(p[0], " / 100%")), bquote(paste(p[c], " / 100%")), bquote(kappa))) + 
  theme(legend.position=c(0.77, 0.9), legend.text.align = 0) +
  theme_td() +
  coord_flip()


