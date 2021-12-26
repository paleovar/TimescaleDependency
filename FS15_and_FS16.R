source("helpers/init.R")

res <- readRDS("data/lr_mle_comp_N200_samples6000.Rds")
beta_seq <- seq(-1, 3, by=0.5)
rmse <- res %>% unnest(data) %>% mutate(dev = ((-1)*beta-slope)**2) %>% group_by(beta, method) %>% summarise(mean=mean(dev), sd_mean=mean(slopesd)) %>% filter(method %in% c("MLE_Rlaw", "LR"))

ggplot(rmse, aes(x=beta)) + geom_point(aes(y=sqrt(mean), shape=method), size=2) + theme_td() +
  scale_x_continuous(limits = c(min(beta_seq), max(beta_seq)), breaks=c(rev(beta_seq)), labels=c(rev(beta_seq)),
    expand = c(0.05, 0.05), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(limits = c(0, 0.30), 
    expand = c(0.01, 0.01), sec.axis = dup_axis(name = NULL, labels = NULL)) +
   theme(panel.background = element_rect(fill = "transparent", colour = "black", size=.8), legend.position = c(0.8, 0.78)) +
  guides(shape=guide_legend(ncol=1)) +
  xlab(TeX("$\\beta$")) + ylab(TeX("root-mean-squared error")) +
  #geom_point(aes(y=rep(rmse$sd_mean[!is.na(rmse$sd_mean)],each=2)), shape=4,
   #             color = "grey50", size=1.5) + 
#  geom_line(aes(x=rep(unique(rmse$beta), each=2),y=rep(rmse$sd_mean[!is.na(rmse$sd_mean)],each=2)), linetype="dashed", color = "grey50", size=0.8) +
  #annotate(geom="text", x=min(beta_seq)+.25, y=mean(rmse$sd_mean, na.rm=T)+0.01, label=TeX("$\\Delta \\beta_{LR}$"),
  #            color="grey30") + 
  theme(legend.background = element_rect(fill="white",
                                  size=0.5, linetype="solid", 
                                  colour ="grey50"),
        legend.title = element_text()) + 
  labs(shape="estimator") +
  scale_shape(solid=T) +
  scale_shape_manual(values=c(17, 15), breaks=c("LR", "MLE_Rlaw"), labels=c("LR", "MLE"))
    
res <- readRDS("data/lr_mle_comp_beta_irr_N200.Rds")

tbb1 <- res %>% unnest(data) %>% group_by(Name, beta) %>% 
  dplyr::summarise(slope_lr_sd=sd(slope_lr),
                   slope_glr_sd=sd(slope_glr),
                   slope_mle_Rlaw_sd=sd(slope_mle_Rlaw)) %>%
  pivot_longer(., cols=c(slope_lr_sd, slope_glr_sd, slope_mle_Rlaw_sd), names_to="method", values_to="beta_est_sd") %>%
  mutate(method=case_when(method=="slope_mle_Rlaw"| method=="slope_mle_Rlaw_sd" ~ "MLE (Rlaw)",
                          method=="slope_lr"| method=="slope_lr_sd" ~ "LR",
                          method=="slope_glr"| method=="slope_glr_sd" ~ "GLR",
                            TRUE ~ method))

tbb2 <- res %>% unnest(data) %>% group_by(Name, beta) %>% 
  dplyr::summarise(slope_lr=mean(slope_lr), 
                   slope_glr=mean(slope_glr),
                   slope_mle_Rlaw=mean(slope_mle_Rlaw),
                   slopesd=mean(slopesd)) %>%
  pivot_longer(., cols=c(slope_lr, slope_glr, slope_mle_Rlaw), names_to="method", values_to="beta_est") %>% 
  mutate(slopesd = ifelse(method != "slope_lr", NA, slopesd)) %>% 
  mutate(method=case_when(method=="slope_mle_Rlaw"| method=="slope_mle_Rlaw_sd" ~ "MLE (Rlaw)",
                          method=="slope_lr"| method=="slope_lr_sd" ~ "LR",
                          method=="slope_glr"| method=="slope_glr_sd" ~ "GLR",
                            TRUE ~ method))

tbb <- inner_join(tbb1, tbb2)

pages.meta <- readRDS("data/pages_meta_scaling.Rds")
ggplot(tbb %>% inner_join(., pages.meta) %>% filter(!method=="GLR") %>% arrange(ID) %>% 
         mutate(method=case_when(method=="MLE (Rlaw)" ~ "MLE", TRUE ~method))
         ) + 
  geom_point(aes(x=beta, y=beta_est, color=as.factor(ID)), size=2, alpha=0.5) + 
  geom_boxplot(aes(x=beta, y=beta_est, group=beta), outlier.shape = NA) + facet_wrap(~method) +
  theme_bw() +
  scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL), breaks=seq(-1, 4, 1), expand=c(0,0), limits=c(-2, 4)) +
  theme_td() +
  xlab(TeX("$\\beta_{sim}$")) + ylab(TeX("$\\beta_{est}$")) + 
  guides(color=guide_legend(ncol=1,byrow=TRUE, name="PAGES ID")) + 
  theme(panel.spacing = unit(.0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.3), 
        legend.text=element_text(size=6),
        legend.spacing.y = unit(0.05, 'cm'),
        panel.grid.major = element_line(colour = "grey")) 
