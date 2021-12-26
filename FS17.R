source("helpers/init.R")
source("helpers/functions.R")
library(nest)

pages.prxlist <- readRDS("data/proxylist.Rds")
tbb <- readRDS("data/beta_N100.Rds")
pages.prxlist <- ltot(pages.prxlist) %>% inner_join(., tbb %>% filter(signal=="pages2k"), by="Name") 

N <- length(pages.prxlist$data)

mat.nexcf_ci <- matrix(data=NA, nrow=N, ncol=N)
mat.nexcf_p <- matrix(data=NA, nrow=N, ncol=N)

for(i in 1:N){
  print("i=")
  print(i)
  for(j in 1:N){
    if(i==j){
      next
    }
    print(j)
    t_i <- tseries::na.remove(pages.prxlist$data[[i]])
    t_j <- tseries::na.remove(pages.prxlist$data[[j]])
    tryCatch({
      mat.nexcf_ci[i,j] <- as.numeric(nexcf_ci(t_i, t_j)$rxy, conflevel=0.05)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    tryCatch({
      mat.nexcf_p[i,j] <- as.numeric(nexcf_ci(t_i, t_j)$pval, conflevel=0.05)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

library(ggcorrplot)

ggcorrplot(mat.nexcf_ci, type = "lower",
           outline.col = "white", pch.cex =1,
           p.mat = mat.nexcf_p, sig.level = 0.05) + 
  theme_td(8) + xlab("PAGES ID") + ylab("PAGES ID") +
  scale_x_continuous(breaks=1:dim(mat.nexcf_ci)[[1]], labels=pages.prxlist$ID,
                     expand = c(0.03, 0.5), 
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks=1:dim(mat.nexcf_ci)[[1]], labels=pages.prxlist$ID,
                     expand = c(0.04, 0.01), 
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low="#4444DE", mid="white", high="#FF7100", 
                       midpoint=0,    
                       breaks=c(-0.45, -0.3, -0.15, 0, 0.15, 0.3, 0.45),
                       limits=c(-0.45,0.45)) +
  theme(legend.position=c(0.15,0.75)) +
  guides(fill = guide_colorbar(barwidth = 0.4, barheight = 5, title="test", position=c(0.1, 0.8)))

print(mean(mat.nexcf_ci, na.rm=T))
print(quantile(unlist(mat.nexcf_ci), na.rm=T,  probs = c(.05, .5, .95)))
