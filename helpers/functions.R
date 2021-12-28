#' @title Reverse log transform
#'
#' @description everse log transformation for plotting in logarithmic space. Useful when plotting power spectral densities.
#' adapted from the metR package https://rdrr.io/cran/metR/src/R/reverselog_trans.R by Elio Campitelli
#' Part of PTBox https://palmodapp.cloud.dkrz.de/index.php/ptbox/
#'
#' @param base Base of the logarithm, default value 10
#'
#' @examples 
#' # 
#' ggplot(t, aes(p, t)) + 
#' geom_line() +
#' coord_flip() +
#' scale_x_continuous(trans = "reverselog")
#' 
#' the base can easily changed: 
#' ggplot(t, aes(p, t)) + 
#' geom_line() +
#' coord_flip() +
#' scale_x_continuous(trans = reverselog_trans(10))
#' #'
#' @family ggplot2 helpers
#' @export
reverselog_trans <- function(base = 10) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
            scales::log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#--------------------------------------------------------------#

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
#' @description plots the power spectral density in log-log space.
#' Part of PTBox https://palmodapp.cloud.dkrz.de/index.php/ptbox/
#' @param tibble nested tibble (name, data) with column Spec(freq, spec, dof, lim.1, lim.2) is the nested spectrum
#' @param main title
#' @param name.y y-axis label
#' @param name.x x-axis label
#' @param name.data name of the nested spectrum column, i.e. either "Spec" or "Smooth_spec"
#' @return ggplot
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
#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
#' @export
GetVar <- function (spec, f, dfreq = NULL, df.log = 0, bw = 3) 
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
#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
#' @export
colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}

#--------------------------------------------------------------#
#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
#' @export
ltot <- function(list_spec_obj){
  
  df <- tibble(Name=character(), data=list())
  for (i in 1:length(list_spec_obj)){ 
    tmp <-  coredata(list_spec_obj)
    if(is.null(names(list_spec_obj))){
      list_spec_obj <- as.list(list_spec_obj)
      names(list_spec_obj) <- as.character(seq(1, length(list_spec_obj),1))}
    df <- df %>% add_row(
      Name = names(list_spec_obj)[[i]], 
      data = list(tmp[[i]])
    )
  }
  return(df)
}

#--------------------------------------------------------------#
#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
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
  ) #%>% select(model, interp.res, EquiTS)
  return(y_equi)
}

#--------------------------------------------------------------#
#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
#' @export
tibble_spec <- function(y, k, nw){
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
                                   restmp <- SpecMTM(target, k=k, nw=nw, detrend=TRUE)
                                   tmpdf[[i]] <- restmp$spec
                                 } 
                                 tmpdf[["freq"]] <- restmp$freq
                                 tmpdf[["dof"]] <- restmp$dof
                                 tmpdf %>% as_tibble()
                               }
                               )
  ) #%>% select(model, interp.res, Spec)
  return(y_spec)
}

#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
#' @export
tibble_to_list <- function(tibble_spec_obj){
  tmp <- as.list(tibble_spec_obj)
  #names(tmp$Spec) <- tmp$model
  tmp$model <- NULL
  tmp$interp.res <- NULL
  lapply(tmp$Spec, function(x){
    x <- as.list(x)
    if ("temp" %in% names(x)){
      x$spec <- x$temp
      x$temp <- NULL
    }
    if ("dR" %in% names(x)){
      x$spec <- x$dR
      x$dR <- NULL
    }
    x$temp <- NULL
    class(x) <- "spec"
    return(x)
  }
  )
}

#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
#' @export
list_to_tibble <- function(list_spec_obj){
  
  df <- tibble(model=character(), data=list()) 
  
  for (i in 1:length(list_spec_obj)){
    
    if(class(list_spec_obj[[i]])=="zoo"){
      time <- index(list_spec_obj[[i]]) #convert to CE
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

#' @title
#' @description 
#' @param
#' @param
#' @param 
#' @return 
#' @export
cut <- function(target, from, to, index=FALSE)
{
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
