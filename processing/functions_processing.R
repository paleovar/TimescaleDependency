#model parameters
names_echam <- c("ECHAM5", "Trace21k", "MPI-M")
echam_weights <- 0.
control_runs <- c("CESM_LM_cont", "MPI-M_cont", "Trace21k_orb")
control_runs_weights <- 0

#--------------------------------------------------------------#
#processing parameters
k=3
nw=2
gregorian_year <- 365.2425
get_min.res <- function(tscales=c(200,10)){min(tscales)}
get_min.range <- function(tscales=c(200,10)){max(tscales)*3}
get_max.hiat <- function(tscales=c(200,10)){min(tscales)*5}

tscale <- list()
tscale$cen <- c(200, 10)
#--------------------------------------------------------------#
#' @title
#' @description 
#' @param target spectrum that needs to be cutted 
#' @param from starting point
#' @param to end point
#' @param index FALSE / TRUE indicates whether "from" and "to" refer to an index or a period (in years)
#' @return object of class "spec"
#' @export
get_restrictions <- function(x){
  if(x == "strong"){
    min.res = 20
    max.hiat = 40 #floor(min.res*2)
    length.min = 30
    min.range= 30 #length.min
  }
  if(x == "loose"){
    min.res = 80
    max.hiat = 160 #floor(min.res*2)
    length.min = 20
    min.range= 20 #length.min
  }
  return(list("min.res"=min.res, "max.hiat"=max.hiat, "min.range"=min.range, "length.min"=length.min))
}

#--------------------------------------------------------------#
#' @title Extract part of spectrum
#' @description Extract part of spectrum
#' @param target spectrum that needs to be cutted 
#' @param from starting point
#' @param to end point
#' @param index FALSE / TRUE indicates whether "from" and "to" refer to an index or a period (in years)
#' @return object of class "spec"
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

#--------------------------------------------------------------#
#' @title get nth element of a vector
#' @description nth element of a vector
#' @param vector target object
#' @param starting_postition integer
#' @param n nth element counted from starting poistion
#' @return numeric() 
#' @export
nth_element <- function(vector, starting_position, n) { 
  vector[seq(starting_position, length(vector), n)] 
}

#--------------------------------------------------------------#
#' @title Spline interpolation of the spectrum
#' @description  Spline interpolation of the spectrum, similar to PaleoSpec::SpecInterpolate()
#' @param freqRef target frequencies (vector)
#' @param spec object of class "spec" with freq, spec and dof list
#' @return object of class spec
#' @export
SpecInterpolateSpline <- function (freqRef, spec) {
  result <- list()
  result$freq <- freqRef
  result$spec <- spline(x=spec$freq, y=spec$spec, xout=freqRef)$y
  result$dof <- spline(x=spec$freq, y=spec$dof, xout=freqRef)$y
  class(result) <- "spec"
  return(result)
}


#--------------------------------------------------------------#
#' @title 
#' @description  
#' @param freqRef
#' @param spec
#' @return
#' @export
transferSpec <- function (specList, input=input.spec, output=output.spec, iRemoveLowest = 1){
  remove.lowestFreq <- function (spec, iRemove){ 
    {
      if (iRemove == 0) 
        index = seq(spec$spec)
      else index <- (-(1:iRemove))
      spec$spec <- spec$spec[index]
      spec$freq <- spec$freq[index]
      spec$dof <- spec$dof[index]
      return(spec)
    }
  }
  get.fend.existing <- function (x){
    return(max(x$freq[!is.na(x$spec)]))
  }
  get.fstart.existing <- function (x) {
    return(min(x$freq[!is.na(x$spec)]))
  }
  get.df <- function (x){
    return(mean(diff(x$freq)))
  }
  AddConfInterval_Fdist <- function(transferspec, var.dof1, var.dof2, pval = 0.05){ 
    {
      if (!(length(transferspec$spec) == length(var.dof1)) && (!length(transferspec$spec) == 
                                                               length(var.dof2))) {
        stop("same lengths must be provided")
      }
      if (!(is.numeric(transferspec$spec)) || !(is.numeric(var.dof1)) || 
          !is.numeric(var.dof2)) {
        stop("non-numeric arguments")
      }
      res <- matrix(NA, nrow = length(transferspec$spec), ncol = 2)
      for (i in 1:length(transferspec$spec)) {
        QF <- qf(p = c(pval/2, (1 - pval/2)), df1 = var.dof1[i], 
                 df2 = var.dof2[i])
        tmp <- QF * transferspec$spec[i]
        res[i, ] <- tmp
      }
      transferspec$lim.1 <- res[,2]
      transferspec$lim.2 <- res[,1]
      class(transferspec) <- "spec"
      return(transferspec)
    }
  }
  specList <- lapply(specList, remove.lowestFreq, iRemove = iRemoveLowest)
  freqRef <- seq(from = min(unlist(lapply(specList, get.fstart.existing))), 
                 to = max(unlist(lapply(specList, get.fend.existing))), 
                 by = min(unlist(lapply(specList, get.df))))
  specList.interpolated <- list()
  for (i in 1:length(specList)) specList.interpolated[[i]] <- SpecInterpolate(freqRef, 
                                                                              specList[[i]])
  NSpectra <- length(specList.interpolated)
  result <- list(freq = specList.interpolated[[1]]$freq, spec = rep(0, 
                                                                    length(specList.interpolated[[1]]$spec)))
  specMatrix <- matrix(NA, NSpectra, length(specList.interpolated[[1]]$spec))
  dofMatrix <- matrix(NA, NSpectra, length(specList.interpolated[[1]]$spec))
  for (i in 1:length(specList.interpolated)) {
    if (sum((result$freq - specList.interpolated[[i]]$freq)^2) > 
        0.1) 
      stop("Different spectra length or resolutions")
    specMatrix[i, ] <- specList.interpolated[[i]]$spec
    dofMatrix[i, ] <- specList.interpolated[[i]]$dof
  }
  
  var.dof1 <- dofMatrix[output, ]
  var.dof2 <- dofMatrix[input,]
  result$spec <- mapply('/', specMatrix[output, ], specMatrix[input, ]) #na.rm=TRUE
  result <- AddConfInterval_Fdist(result, dofMatrix[output, ], dofMatrix[input,])
  class(result) <- "spec"
  return(list(spec = result))
}
