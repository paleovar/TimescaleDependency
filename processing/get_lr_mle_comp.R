source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
library(foreach)
library(doParallel)
library(furrr)
library(poweRlaw)

RlawFit <- function(x){
  m <- conpl$new(x)
  m$setXmin(estimate_xmin(m))
  return(-1/(estimate_pars(m)$pars-1))
}

scale <- "cen"
ncores <- detectCores()

min.res = get_min.res(tscale[[scale]])
min.range = get_min.range(tscale[[scale]])
max.hiat = get_max.hiat(tscale[[scale]])

save <- F

#### Case 1: Regular data ("data/supp/lr_mle_comp_N=",N, "_samples", samples,".Rds")
#Comparison of linear regression and maximum likelihood estimation of scaling exponents for irregularly sampled data**
#The following methods are compared:
#- **Linear regression (LR)** from the `PaleoSpec::SlopeFit()` function based on the stats::lm() function after log-binning. 
#- **Maximum likelihood estimation (MLE)** as suggested by Clauset, Shalizi, Newman (2007). We use both the `poweRlaw::estimate_pars()` (denoted by "MLE_Rlaw").

#We simulate surrogate time series with power-law scaling $\beta$ using `PaleoSpec::SimPowerlaw()`. 
#Afterwards, we irregularly sample them using a constant mean resolution with fixed variance. 
#We apply block-averaging to the new resolution, compute the spectrum using the multitaper method `PaleoSpec::SpecMTM()` 
#and estimate the scaling coefficient again. 
#The statistics of the two methods **LR** and **MLE** are compared to assess their performance.

beta_lr_mat <- tibble(
  beta = numeric(),
  slope = numeric(),
  slopesd = numeric()
)

res_lr <- tibble()
res_mle_Rlaw <- tibble()

beta_mle_Rlaw_mat <- beta_lr_mat

N <- 200 #! Caution with these numbers. Compuing time and memory demand might be high
samples <- 6000 #! Caution with these numbers. Compuing time and memory demand might be high
#mean sample lengths of proxies
#lengths <- tbb[[scale]] %>% mutate(len = purrr::map(data, function(.data){length(.data)})) %>% select(len)
#samples <- mean(unlist(lengths)) #(this is 900)

### Mean-squared error of estimator
#### Case 1: Regular data
#We compute the performance of both estimators and evaluate it using the measure of mean-squared error  ${\displaystyle \operatorname {MSE} ={\frac {1}{n}}\sum _{i=1}^{n}(Y_{i}-{\hat {Y_{i}}})^{2}.}$.

#**1.1. Generate data**
myCluster <- makeCluster(ncores, # number of cores to use
                         type = "PSOCK")
registerDoParallel(myCluster)

beta_seq <- seq(-1, 3, by=0.5)

for(beta in beta_seq){
  print(beta)
  powerlaw <- foreach(i = 1:N) %dopar% {
    as.ts(PaleoSpec::SimPowerlaw(beta, samples))
  }
  
  spec <- mclapply(powerlaw, function(x) SpecMTM(x, detrend=T), mc.cores=ncores)
    
  #linear regression estimates
  beta_lr <- mclapply(spec, function(x) PaleoSpec::SlopeFit(x, 1/(samples/2), 1/2, bDebug = F), mc.cores=ncores)
  
  
  #maximum likelihood estimates
  beta_mle_Rlaw <- mclapply(spec, function(x) RlawFit(x$spec), mc.cores=ncores)
  
  #summary linear regression
  res_lr <- rbind(res_lr, beta_lr_mat %>% 
                    add_row("beta" = beta, "slope" = unlist(lapply(beta_lr, function(x) x$slope)), 
                          "slopesd" = unlist(lapply(beta_lr, function(x) x$slopesd))) %>% 
                    group_by(beta) %>% nest() %>% ungroup())
  
  #summary maximum likelihood
  res_mle_Rlaw <- rbind(res_mle_Rlaw, beta_mle_Rlaw_mat %>% add_row("beta" = beta, "slope" = unlist(beta_mle_Rlaw)) %>% group_by(beta) %>% nest() %>% ungroup())
  
}

stopCluster(myCluster)

res <- rbind(res_lr %>% add_column("method"= "LR"), res_mle_Rlaw %>% add_column("method"= "MLE_Rlaw"))

if(save){
  saveRDS(res, ("data/supp/lr_mle_comp_N=",N, "_samples", samples,".Rds"))
}

#### Case 2: Irregular data ("data/supp/lr_mle_comp_beta_irr_N",N, ".Rds")
#**2.1. Generate data**
prxlist <- readRDS("data/proxylist.Rds")
df <- prxlist %>% mutate(diffs=purrr::map(data, function(x) diff(index(x)))) %>% select(-data)
specs <- readRDS("data/proxy_spectra.Rds")
tbb_init <- inner_join(tibble(Name= specs[["Name"]]), 
                specs %>% select(Name, interp.res)) %>% inner_join(., df)

betas <- c(-1,-0.5, 0.,0.5, 1,1.5, 2,2.5, 3) #testing the code with a few values for beta is recommended
N <- 200 #!Caution with this number, high memory and computing power might be needed

myCluster <- makeCluster(ncores, # number of cores to use
                         type = "PSOCK")
registerDoParallel(myCluster)

start.time <- Sys.time()
t_final <- foreach(j = betas, .packages=c('foreach', 'tidyr', 'dplyr', 'tibble')) %dopar% {
  beta <- j
  print(beta)
  result <- foreach(i = 1:N, .packages=c('dplyr', 'tibble', 'purrr', 'zoo', 'PaleoSpec', 'tseries', 'RScaling', "poweRlaw", 'furrr')) %dopar% {
      if(i %in% seq(1, 1000, by=100)){print(i)}
      #create artificial proxy using block-averaging to sample solution
      tbb <- tibble()
      tbb <- tbb_init %>% add_column("beta"=beta) %>% 
        mutate(sample = future_map2(diffs, beta, function(.diffs, .beta) sample_from_proxy(.diffs, .beta, plt=F)[["sample"]])) 
      
      #convert to time series using equidistant
      tbb <- tbb %>% mutate(ts = future_map2(sample, interp.res, function(.sample, .interp.res) 
                         MakeEquidistant(index(.sample), coredata(.sample), dt=.interp.res))) 
      
      #compute spectrum
      tbb <- tbb %>% mutate(spec = future_map(ts, function(.ts){SpecMTM(na.remove(.ts), k=k, nw=nw, detrend=TRUE)})) 
      
      #fit scaling coefficients to spectrum (spec)
      tbb <- tbb %>% 
                      #lr method
                      mutate(temp = future_map(spec, function(.spec){SlopeFit(.spec, min(.spec$freq)*2, max(.spec$freq)/2, bDebug=F)})) %>%
                                         mutate(slope_lr = future_map(temp, function(.temp){.temp$slope})) %>% 
                                         mutate(slopesd = future_map(temp, function(.temp){.temp$slopesd})) %>%
                                         select(-temp) %>% 
                      #mle_Rlaw method
                      mutate(slope_mle_Rlaw = future_map(spec, function(.spec){RlawFit(.spec$spec)})) 
       
      tibble(Name=tbb$Name, 
             beta=tbb$beta,
             slope_lr=-1*unlist(tbb$slope_lr), 
             slope_mle_Rlaw=-1*unlist(tbb$slope_mle_Rlaw),
             slopesd = unlist(tbb$slopesd),
             )
    }
    results <- tibble()
    for(i in 1:length(result)){results <- rbind(results, result[[i]])}
    results
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
stopCluster(myCluster)

print("start rbind")
res <- tibble()
for(i in 1:length(t_final)){res <- rbind(res, t_final[[i]])}

print("start nest")
res <- res %>% group_by(Name, beta) %>% tidyr::nest() %>% ungroup()

if(save){
  saveRDS(res,  paste0("data/supp/lr_mle_comp_beta_irr_N",N, ".Rds"))
}
