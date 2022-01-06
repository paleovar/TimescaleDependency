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
