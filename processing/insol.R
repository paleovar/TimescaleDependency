library(palinsol)
library(rlist)
#time in years after 1950

#standard insolation function from palinsol package
insolation <- function(times, astrosol=ber78, long=pi/2, lat=65*pi/180){
  sapply(times, function(tt) Insol(orbit=astrosol(tt), long=long, lat=lat))
  }

#adapted i.e. taken from https://bitbucket.org/mcrucifix/insol
Insol_H <- function(orbit,long, lat=65*pi/180,S0=1365, H=NULL)
{
  # Returns daily mean or hourly incoming solar insolation after Berger (1978) as function of the hour angle of the sun
  
  # eps  : obliquity (radians)
  # varpi: longitude of perihelion (radians)
  # e : eccentricity 
  # long : true solar longitude 
  #        (radians; = pi/2 for summer solstice )
  #                     pi   for autumn equinox  )
  #                     3* pi/2 for winter solstice )
  #                     0 for spring equinox )
  # lat : latitude
  # orbit is a list containing eps, varpi and e (see Ber78)
  # S0 : Total solar irradiance (default : 1365 W/m^2)
  # returns : daily mean incoming solar radiation at the top of the atmosphere
  #           in the same unit as S0. 
  # H : if NULL (default): compute daily mean insolation
  # else: hour angle (in radians) at which insolation is being computed
  
  varpi <- NULL
  eps <- NULL
  ecc <- NULL
  for (i in names(orbit)) assign(i,orbit[[i]])
    nu <- long - varpi
    rho <- (1-ecc^2)/(1+ecc*cos(nu))
    sindelta <- sin(eps)*sin(long)
    cosdelta <- sqrt(1-sindelta^2)
    sinlatsindelta <- sin(lat)*sindelta
    coslatcosdelta <- cos(lat)*cosdelta 
  if (is.null(H))
  {
    cosH0 <- min(max(-1,-sinlatsindelta/coslatcosdelta),1)
    sinH0 <- sqrt(1-cosH0^2)
    H0 <- acos(cosH0)
    print(H0)
    insol <- S0/(pi*rho^2)*(H0*sinlatsindelta+coslatcosdelta*sinH0)
  }
  else 
  {
    insol <- pmax(0, S0/(rho^2)*(sinlatsindelta+coslatcosdelta*cos(H)))
  }
  return(insol)
}

insolation_year <- function(times, astrosol=ber78, long=pi/2, lat=65*pi/180){sapply(times, function(tt) Insol_H(orbit=astrosol(tt), long=long, lat=lat,  H=NULL))}

insolation_day <- function(times, astrosol=ber78, lat=65*pi/180)
{ 
  seqs <- seq(from=120, to=360, by=120)
  dailyTSI <- lapply(times, function(tt){Insol(astrosol(tt), long = day2l(astrosol(tt), seqs), lat = lat)})
  dailyTSI.mat <-list.rbind(dailyTSI)
  vec <- as.vector(t(dailyTSI.mat))
  dat <-zoo(vec, order.by = seq(from = times[1], to = times[length(times)]+1, by = 120/360)[-1])
  return(dat)
}

insolation_hour <- function(times, astrosol=ber78, tmin=0.5, lat=65*pi/180)
{
  hourly_insolation <- function(tt, astrosol=ber78, tmin, lat = 65*pi/180){
    hours <- seq(from=-12, to=12-tmin, by=tmin)
    hourlyTSI <- lapply(seq(from = 1, to = 360, by = 1), function(d){
      Insol_H(astrosol(tt), long = day2l(astrosol(tt), d), lat = lat, H=hours*pi/12)
      })
    hourlyTSI.mat <-list.rbind(hourlyTSI)
    vec <- as.vector(t(hourlyTSI.mat))
    dat <-zoo(vec, order.by =seq(from = tt, to = tt+1, by = 1/360*(diff(hours)[1]/24))[-1])
    return(dat)}
  years_list <- lapply(times, function(tt) hourly_insolation(tt, astrosol, tmin=tmin, lat=lat))
  list.rbind(years_list)
}

SpecInterpolateSpline <- function(freqRef, spec) 
{
  result <- list()
  result$freq <- freqRef
  result$spec <- spline(x=spec$freq, y=spec$spec, xout=freqRef)$y
  result$dof <- spline(x=spec$freq, y=spec$dof, xout=freqRef)$y
  class(result) <- "spec"
  return(result)
}

#days
INS <- list()
INS$tts <- seq(from = -10000, to = 50, by = 1) #take care of time units
INS$summer<- insolation_day(times = INS$tts, astrosol= ber78, lat=65*pi/180) 
x_d <- as.ts(INS$summer)
smooth.spec <- PaleoSpec::SpecMTM(x_d, k=3, nw=2)
spec.df <- as.data.frame(smooth.spec$freq, spec = smooth.spec$spec)
l <- list()
l$log_year <- LogSmooth(smooth.spec, df.log = 0.001, removeFirst = 1, removeLast = 30)
l$spec_year <- smooth.spec
#LPlot(cut(l$log_year, 1000, 0.75, index=FALSE))
l$log_year_cutted <- cut(l$log_year, 1000, 0.75, index=FALSE)


#hours
INS <- list()
INS$tts <- seq(from = 20, to = 40, by = 1) #take care of time units
INS$summer<- insolation_hour(times = INS$tts, astrosol= ber78, tmin=4, lat=65*pi/180) 
x_d <- as.ts(INS$summer)
smooth.spec <- PaleoSpec::SpecMTM(x_d, k=3, nw=2)
l$log_hours <- LogSmooth(smooth.spec, df.log = 0.002, removeLast = 1000, removeFirst = 1)
l$spec_year <- smooth.spec
#LPlot(l$log_hours)
#LPlot(cut(l$log_hours, 1.8/360, 1/(360*1.5), index=FALSE))
l$log_hours_cutted <- cut(l$log_hours, 1.8/360, 1/(360*1.5), index=FALSE)

s <- list()
s$hours <- l$log_hours_cutted
s$years <- l$log_year_cutted
  
res <- MeanSpec(s)
result <- res$spec

#if necessary, spectra can be further interpolated using PaleoSpec::SpecInterpolate() or SpecInterpolateSpline()
