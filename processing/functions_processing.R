##cut_warming <- function(tibble, cut_time=2020, cut=F){
#  if(cut==FALSE){return(tibble)}
#  if(cut==TRUE){
#    tibble <-  rowid_to_column(tibble, "ID")
#    tbb1 <- tibble %>% filter(signal %in% c(signal_tbb %>% filter(type=="obs"&signal!="pages2k") %>% select(signal))$signal|grepl("highres",signal))
#    tbb2 <- tibble %>% filter(!signal %in% tbb1$signal) %>% 
#      mutate(data = purrr::map(data, ~ filter(., time <= cut_time)))
#    res <- rbind(tbb1, tbb2) %>% arrange(ID) %>% select(-ID)
#    return(res)
#  }
#}

#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
#' @export
cut <- function(target, from, to, index=FALSE){
  if (index==FALSE){
   target <- lapply(target, function(x) x[max(which.min(abs(1/target$freq - from)), 0):min(which.min(abs(1/target$freq - to)), length(target$freq))])
    class(target) <- "spec"
    return(target)
  }
  if (index==TRUE)
  {
    target <- lapply(target, function(x) x[from:to])
    class(target) <- "spec"
    return(target)
  }
}


nth_element <- function(vector, starting_position, n) { 
  vector[seq(starting_position, length(vector), n)] 
}

SpecInterpolateSpline <- function (freqRef, spec) {
  result <- list()
  result$freq <- freqRef
  result$spec <- spline(x=spec$freq, y=spec$spec, xout=freqRef)$y
  result$dof <- spline(x=spec$freq, y=spec$dof, xout=freqRef)$y
  class(result) <- "spec"
  return(result)
}
