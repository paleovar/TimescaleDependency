source("helpers/init.R")
source("helpers/functions.R")
source("processing/functions_processing.R")
speclist_smoothed_tbb <- readRDS("data/global_mean_spectra.Rds")

obs <- c(signal_tbb %>% filter(type=="obs"))$signal

w <- rep(1, length(speclist_smoothed_tbb$signal))
w[which(speclist_smoothed_tbb$signal %in% names_echam)] <- echam_weights
w[which(speclist_smoothed_tbb$signal %in% control_runs)] <- control_runs_weights
w[which(speclist_smoothed_tbb$signal %in% obs)] <- 0


speclist_smoothed <- tibble_to_list(speclist_smoothed_tbb, name.data="data")
names(speclist_smoothed) <- speclist_smoothed_tbb$signal

speclist_smoothed[which(w==0)] <- NULL
names(speclist_smoothed)
speclist_smoothed_tbb$mean <- NULL

w <- rep(1, length(names(speclist_smoothed)))

speclist <- lapply(speclist_smoothed, function(x) cut(x, from=5000, 12/9, index=FALSE))
l <- list()
len <- length(w[w!=0])
for(i in 1:1000){
  w <- sample(0:1000, len, replace = TRUE)/1000
  l[[i]] <- MeanSpec(speclist, weights=w)$spec 
  #l$mean$dof <- NULL
  l[[i]]$lim.1 <- NULL
  l[[i]]$lim.2 <- NULL
}

LPlot(l[[20]])
length(l)
