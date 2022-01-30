if(!exists("samples")){
  samples <- list()
}


# LARGELY TAKEN FROM: Neukom, Raphael; Barboza, Luis A.; Erb, Michael P.; Shi, Feng; Emile-geay, Julien; Evans, Michael N.; et al. (2019): 
# Global mean temperature reconstructions over the Common Era. figshare. Collection. https://doi.org/10.6084/m9.figshare.c.4507043.v2 
#-------------------------------#

### The original code serves to generate Figures for GMST reconstruction paper: PAGES2k Consortium (2019); Nature Geoscience. RN 2019/03/29

#functions & definitions ------------------------------------------------
### read a table containing a time series, first column = years
read.ts<-function(filename,sep=";",header=T,skip=0){
  ind<-as.matrix(read.table(filename,sep=sep,header=header,skip=skip))
  data<-ts(ind[,-1],start=ind[1,1])
  return(data)
}

### apply a function to a timeseries-matrix. the outcome is a time series with the same tsp
tsapply<-function(x,dim,fun){
  out<-ts(apply(x,dim,fun),start=start(x)[1])
}

#folder with reconstruction outputs
recons.folder <- "processing/raw_data/pages2k/"

experiment.names <-c("CPS","PCR","M08","PAI","OIE","BHM","DA")

#reference period
start.ref<-1770
end.ref<-1850

##read in the reconstructions------------------------
nexp<-length(experiment.names)
recon.files<-paste0(experiment.names,".txt")

recons<-list()
for(i in seq_along(recon.files)){
  filename<-paste0(recons.folder,recon.files[i])
  recons[[i]]<-read.ts(filename,sep="\t",header=F,skip=1)  
}

for(exp in seq_along(recon.files)){
  em<-tsapply(recons[[exp]],1,function(x) median(x,na.rm=T))
  em.ref<-mean(window(em,start.ref,end.ref))
  recons[[exp]]<-recons[[exp]]-em.ref
}

s <- function(x){PaleoSpec::SpecMTM(x, k=3, nw=2, detrend=TRUE)}

specs <- lapply(recons, function(x) apply(x, 2, s))

samples$recons <- specs

rm(specs, recons, recons.folder, experiment.names, em.ref, em, filename)
gc()