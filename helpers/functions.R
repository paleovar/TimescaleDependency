
#' @title Mean Spectrum
#' @description Wrapper for the PaleoSpec::MeanSpectrum function. We introduced a correction to the number of records, 
#' in case spectra where only partially overlapping.
#' Part of PTBox https://palmodapp.cloud.dkrz.de/index.php/ptbox/
#' @param speclist list of spectra
#' @param iRemoveLowest number of lowest frequencies to remove (e.g. to remove detrending bias)
#' @param weights vector of weights (same length as elements in speclist)
#' @return list(spec,nRecords) spec=average spectrum, nRecords = number of records contributing to each spectral estimate
#' @export
MeanSpec <- function(specList, iRemoveLowest = 1, weights = rep(1, length(specList))){
  meanspec <- PaleoSpec::MeanSpectrum(specList, iRemoveLowest, weights)
  meanspec$spec$spec <- meanspec$spec$spec * length(specList)/meanspec$nRecord
  meanspec$spec <- AddConfInterval(meanspec$spec)
  return(meanspec)
}

#--------------------------------------------------------------#

#' @title PSD plot
#' @description Plots the power spectral density in log-log space.
#' Part of PTBox https://palmodapp.cloud.dkrz.de/index.php/ptbox/
#' @param tibble nested tibble (name, data) with column Spec(freq, spec, dof, lim.1, lim.2) is the nested spectrum
#' @param main title
#' @param name.y y-axis label
#' @param name.x x-axis label
#' @param name.data string, name of the nested spectrum column, i.e. either "Spec" or "Smooth_spec"
#' @param name.col string, grouping feature for colored lines
#' @param name.fill string, grouping feature for filled confidence bands
#' @param name.line string, grouping feature for linetype
#' @param name.alpha string, grouping feature for transparence of filled ribbons
#' @param conf.band TRUE / FALSE denotes whether confidence bands of spectra are plotted or not
#' @return ggplot object
#' @export
plot_spec <- function(tibble, ylims, xlims, name.y=TeX('PSD $S(\\tau)\\, (K^2 yr)$ '), name.x=TeX('period $\\tau\\,(yr)$'), name.data = "data", name.col="signal", name.fill=NULL, name.line=NULL, name.alpha=NULL, conf.bands=T){
    if(is.null(name.fill)){name.fill<-name.col}
    xlims <- rev(sort(xlims))
    ylims <- sort(ylims)
    yrs.period <- rev(c(0.0001, 0.001, 0.01,  0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
    yrs.labels <- rev(c(TeX('$10^{-4}$'),TeX('$10^{-3}$'), TeX('$10^{-2}$'), TeX('$10^{-1}$'), TeX('$10^{0}$'), TeX('$10^{1}$'), TeX('$10^{2}$'), TeX('$10^{3}$'), TeX('$10^{4}$'), TeX('$10^{5}$'), TeX('$10^{6}$')))
    tmp <- tibble %>% unnest(name.data) %>%
        ggplot(aes(x = 1/freq, y = spec)) 
    if(conf.bands==T){
      if(!is.null(name.alpha)){
      tmp <- tmp + geom_ribbon(aes_string(ymin="lim.1", ymax="lim.2", fill=name.fill, alpha=name.alpha), linetype = 0) + 
      guides(fill=FALSE)
    } else {
      tmp <- tmp + geom_ribbon(alpha=0.2, aes_string(ymin="lim.1", ymax="lim.2", fill=name.fill), linetype = 0) + 
    guides(fill=FALSE)  
    }
    }
    tmp + geom_line(size=0.4, aes_string(color=name.col, linetype=name.line)) + theme_td() +
        scale_y_log10(name=name.y, label = scales::trans_format("log10", scales::math_format(10^.x)), expand=c(0.0, 0.0), limits=ylims, sec.axis = dup_axis(name = NULL, labels = NULL))  +
        scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, name=name.x,expand=c(0.0, 0.0), limits=xlims, sec.axis = dup_axis(name = NULL, labels = NULL)) 
    }

#--------------------------------------------------------------#
#' @title Variance computation by integration of the spectrum
#' @description Simple computation of variance on frequency intervals by integration of the spectrum, akin to PaleoSpec::GetVarRatio(), 
#' but without consideration of degrees of freedom not needed here
#' @param spec tibble or list with freq and spec column
#' @param f target frequency interval c(f1,f2) 
#' @param dfreq default NULL, if not NULL spectra are interpolated with dfreq before computing the variance
#' @return  returns variance
#' @export
GetVar <- function (spec, f, dfreq = NULL) 
{
  if (f[1] >= f[2]) 
    stop("f1 must be < f2")
  freqVector <- spec$freq
  if (f[1] < FirstElement(freqVector)) {
    warning("f[1] is smaller than the lowest frequency in the spectrum, set to the lowest frequency")
    f[1] <- FirstElement(freqVector)
  }
  if (f[2] > LastElement(freqVector)) {
    warning("f[2] is larger than the highest frequency in the spectrum, set to the highest frequency")
    f[2] <- LastElement(freqVector)
  }
  Intp <- function (freqRef, spec) 
  {
    result <- list()
    result$freq <- freqRef
    result$spec <- approx(spec$freq, spec$spec, freqRef)$y
    class(result) <- "spec"
    return(result)
  }
  if (is.null(dfreq)) 
    dfreq <- min(diff(spec$freq)[1]/5, (f[2] - f[1])/100)
  newFreq <- seq(from = f[1], to = f[2], by = dfreq)
  vars <- mean(Intp(newFreq, spec)$spec)
  return(vars)
}

#--------------------------------------------------------------#
#' @title Color ramp palette with transparency 
#' @description Color ramp palette with transparency. Taken from https://rdrr.io/github/JBrenn/Helper4me/man/colorRampAlpha.html
#' @param ... arguments to pass to colorRampPalette
#' @param n  	number of colors in palette
#' @param alpha alpha channel (opacity) value [0;1]
#' @return vector of defined colors with transparency interpolating the given sequence
#' @export
colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}


#--------------------------------------------------------------#
#' @title Equidistant timeseries from tibble
#' @description  Wrapper for applying PaleoSpec::MakeEquidistant() to a tibble with data column that contains the timeseries
#' @param y tibble object with nested data column
#' @return tibble object with two nested tibble. The "data" column the raw data and the "EquiTS" contains the equidistant timeseries. 
#' @export
equidistant <- function(y){
  cnt <- 0
  y_equi <- y %>% add_column(EquiTS =
                               .$data %>% lapply(., function(x){
                                 cnt <<- cnt+1
                                 tmp <- x %>% as_tibble() 
                                 tmpdf <- list() 
                                 for (i in names(tmp)){
                                   if(i=="time"|i=="year"){next}
                                   tmpdf[[i]] <- MakeEquidistant(tmp$time, tmp[[i]], dt = y$interp.res[cnt])
                                 }
                                 tmpdf %>% as_tibble()
                               }
                               )
  ) 
  return(y_equi)
}

#--------------------------------------------------------------#
#' @title Computing the spectrum over a tibble
#' @description wrapper for cmputing the spectrum from a tibble using the PaleoSpec pacakge
#' @param y tibble with data EquiTS column that contains interpolated time series
#' @param k a positive integer, the number of tapers, often 2*nw.
#' @param nw a positive double precision number, the time-bandwidth parameter.
#' @return tibble with nested spectrum
#' @export
tibble_spec <- function(y, k=3, nw=2){
  cnt <- 0
  y_spec <- y %>% add_column(Spec =
                               .$EquiTS %>% lapply(., function(x){
                                 cnt <<- cnt+1
                                 tmp <- x %>% as_tibble() 
                                 tmpdf <- list() 
                                 for (i in names(tmp)){
                                   target <- tmp[[i]]
                                   if(any(is.na(tmp[[i]]))){
                                     if("model" %in% names(y)){
                                       print(paste0(y$model[[cnt]], " ", i, " contains NA values, which were removed"))
                                       target = na.remove(tmp[[i]])
                                     }
                                     if("Name" %in% names(y)){
                                       print(paste0(y$Name[[cnt]], " ", i, " contains NA values, which were removed"))
                                       target = na.remove(tmp[[i]])
                                     }
                                   }
                                   target <- target - mean(target)
                                   restmp <- SpecMTM(target, k, nw, detrend=TRUE)
                                   tmpdf[[i]] <- restmp$spec
                                 } 
                                 tmpdf[["freq"]] <- restmp$freq
                                 tmpdf[["dof"]] <- restmp$dof
                                 tmpdf %>% as_tibble()
                               }
                               )
  )
  return(y_spec)
}

#--------------------------------------------------------------#
#' @title Reverse log transform
#' @description reverse log transformation for plotting in logarithmic space. Useful when plotting power spectral densities.
#' adapted from the metR package https://rdrr.io/cran/metR/src/R/reverselog_trans.R by Elio Campitelli
#' Part of PTBox https://palmodapp.cloud.dkrz.de/index.php/ptbox/
#' @param base Base of the logarithm, default value 10
#' @return axis object
#' @export
reverselog_trans <- function(base = 10) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
            scales::log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#--------------------------------------------------------------#
#' @title Tibble spec to speclist converter
#' @description Converts a tibble with nested tibbles that contain the spectrum to a list with objects of class "spec" 
#' @param tibble_spec_obj tibble object that contains the spectrum
#' @param name.data name of the tibble column that contains the spectral information
#' @return list of objects of class "spec"
#' @export
tibble_to_list <- function(tibble_spec_obj, name.data="Spec"){
  tmp <- as.list(tibble_spec_obj)
  tmp$model <- NULL
  tmp$interp.res <- NULL
  res <- lapply(tmp[[name.data]], function(x){
  x <- as.list(x)
  n <- names(x)
  idx <- which(!names(x) %in% c("dof", "freq", "f", "lim.1", "lim.2"))
  n[[idx]] <- "spec"
  names(x) <- n
  class(x) <- "spec"
  return(x)
   }
  )
  return(res)
}

#--------------------------------------------------------------#
#' @title Speclist to tibble spec converter
#' @description Converts a list with objects of class "spec" to tibble with nested tibbles that contain the spectrum. 
#' Inverts the tibble_to_list() functions.
#' @param list_spec_obj  list of objects of class "spec"
#' @return tibble object that contains the spectrum
#' @export
list_to_tibble <- function(list_spec_obj){
  
  df <- tibble(model=character(), data=list()) 
  
  for (i in 1:length(list_spec_obj)){
    
    if(class(list_spec_obj[[i]])=="zoo"){
      time <- index(list_spec_obj[[i]])
      temp <- coredata(list_spec_obj[[i]])
      list_spec_obj[[i]] <- list(time=time, temp=temp)
    }
    
    class(list_spec_obj[[i]]) <- "list"
    
    df <- df %>% add_row(
      model = names(list_spec_obj)[[i]], 
      data = list_spec_obj[[i]] %>% as.data.frame() %>% as_tibble() %>% list()
    )
  }
  return(df)
}

#--------------------------------------------------------------#
#' @title Surrogate Proxies with powerlaw scaling
#' @description Generate proxy surrogates with the same temporal resolution and beta scaling using PaleoSpec::SimPowerlaw() and PaleoSpec::AvgToBin()
#' @param diffs temporal resolution of proxies, e.g. computed from diff(index(zoo object))
#' @param beta scaling coefficient
#' @return list() object with the sample timeseries and simulated powerlaw spectrum
#' @export
sample_from_proxy <- function(diffs, beta){
  #get temporal resolution of proxy records
  #simulate time series with powerlaw scaling, resolution of 1 year and same length as proxy 
  powerlaw <- PaleoSpec::SimPowerlaw(beta, ceiling(sum(diffs)))

  sums <- cumsum(diffs)
  
  avg <- PaleoSpec::AvgToBin(index(powerlaw), powerlaw, breaks= c(sums))
  df <- zoo(avg$avg, order.by = avg$breaks)
  df <- df[1:(length(df)-1)]
  powerlaw <- powerlaw[1:length(powerlaw)-1]
  
  plot(powerlaw, col="red", type="p")
  lines(df, col="blue", type="p")
  
  if(any(unique(index(df)) != index(df))){print("non-unique")}
  return(list(sample=df, powerlaw=powerlaw))
}


##' @title
##' @description 
##' @param
##' @param
##' @param 
##' @return 
##' @export
#cut <- function(target, from, to, index=FALSE)
#{
#  if (index==FALSE){
#    target <- lapply(target, function(x) x[max(which.min(abs(1/target$freq - from)), 0):min(which.min(abs(1/target$freq - to)), length(target$freq))])
#    class(target) <- "spec"
#    return(target)
#  }
#  if (index==TRUE)
#  {
#    target <- lapply(target, function(x) x[from:to])
#    class(target) <- "spec"
#    return(target)
#  }
#}
