######################################################################## 
# frog_power script simulates focalpairs data sets with variable calling 
# activity, sampling randomly from empirical distributions of calling duration
# for each species, and the ratio of green frogs to bullfrogs. calculates 
# p-values for 

#integrates with
# SGE as array job, splitting the job by the number of calls (activity)
############################################################################
rm(list=ls())
library(parallel)

getdata <- function(file) {
  
  fd <- read.csv(file) # import
  
  # convert from two column sets to one column set
  gf=fd[,c("file","onset","duration")]
  bf=fd[,c("file.1","onset.1","duration.1")];
  colnames(bf)=c("file","onset","duration")
  fd=rbind(gf,bf)
  fd=fd[!(fd$duration)=="",]
  
  # Organize and add a few variables
  id=factor(substr(fd$file,2,8))
  species=factor(substr(fd$file,1,1),label=c("B","G"))
  month=as.numeric(substr(fd$file,3,3))
  day=as.numeric(substr(fd$file,4,5))
  minute=as.numeric(substr(fd$onset,1,2))
  second=as.numeric(substr(fd$onset,4,7))
  duration <- as.numeric(substr(fd$duration,4,7))
  
  # calculate onset and offset
  onset=minute*60+second
  offset=onset+duration
  
  fd <- data.frame(id, species, onset, offset, duration)
  fd <- fd[!(duration==0),] # remove 0 durations
  rownames(fd) <- 1:nrow(fd)
  return(fd)
}
getinterrupts <- function(td,random=F){
  gint <- function(l){
    sp <- td$species[l]
    tonset <- td$onset[l]
    # find interrupted calls  (including self)
    inter <- (tonset>=td$onset) & (tonset<=td$offset)  
    inter[l] <- FALSE  #set self to false
    inter <- which(inter)[1]  # get index of first interrupted call
    
    ## get most recent but ended offset time
    #find all previously ending calls
    prev <- which(td$offset<=tonset)  
    # test for conspecific previously ending calls
    hom <- td$species[prev]==td$species[l] 
    # find time between the one most recent previously ending call and current 
    # call (l) for both conspecifics (ptimehom) and heterospecifics (ptimehet)
    pthom <- tonset-td$offset[prev][hom]
    ptimehom <- pthom[which.min(tonset-td$offset[prev][hom])]
    pthet <- tonset-td$offset[prev][!hom]
    ptimehet <- pthom[which.min(tonset-td$offset[prev][!hom])]
    
    if(length(inter)>1) print("interrupted more than one call")
    ttd=data.frame(td[l,],
                   # get species of interrupted call
                   int_sp=ifelse(is.na(inter),"none", 
                                 as.character(td$species[inter])),
                   # get time until interruption
                   itime=ifelse(is.na(inter), NA, tonset-td$onset[inter]),
                   # get interrupt duration
                   int_dur=ifelse(is.na(inter), NA, 
                                  min(td$offset[inter],td$offset[l])-tonset),
                   ptimehet=ifelse(length(ptimehet)>0,ptimehet,NA),
                   ptimehom=ifelse(length(ptimehom)>0,ptimehom,NA))
    return(ttd)
  }
  td2 <- do.call(rbind.data.frame, mclapply(1:nrow(td), gint, mc.cores=8))
  return(td2)}
p_gints <- function(td) {
  n_gcalls <- sum(td$species=="G")
  n_gints <- sum(td$species=="G" & td$int_sp=="B")
  p_gints <- n_gints/n_gcalls
  return(p_gints)
}
ht.stats <- function(data, emp=0, alpha=0.05, enn=NULL,
                     one.sided=TRUE, side="lower") {
  pvalue <- ecdf(data)
  if (one.sided) {
    if (side=="lower") cis <- quantile(data, c(alpha, 1))
    else cis <- quantile(data, c(0, 1 - alpha))
    pv <- ifelse(side=="lower", pvalue(emp), 1 - pvalue(emp))
  } 
  else {
    cis <- quantile(data, c(alpha / 2, 1 - (alpha / 2)))
    pv <- ifelse(emp > mean(data), pvalue(emp)*2, (1 - pvalue(emp)) * 2 )
  }
  return(list(ci=cis, pv=pv))
}


# function computes frequency distrubutions for interruptions at variety
# of calling activity levels (from 50 to 600 calls in an hour). Next, consider 
# increasing number of simulations from 10 to 1000, and making an array job 
# on SGE

fpow <- function(data, ncalls, sims=10) {
  do.call(c, lapply(1:sims, function(x) {
    bdur <- data[data$species=="B", "duration"] # get all bf durations
    gdur <- data[data$species=="G", "duration"] # get all gf durations
    gprop <- do.call(c, lapply(split(data, data$id), function(x) {
      sum(x$species=="G")/nrow(x) # get proportion gf for each file
    }))
    onset.pts <- seq(0, 1798.5, by=.1) 
    gf.calls <- round(sample(gprop, size=1)*ncalls) # number of gf calls
    bf.calls <- ncalls - gf.calls # number of bf calls
    bf.onsets <- sample(onset.pts, size=bf.calls) # draw onset times
    gf.onsets <- sample(onset.pts, size=gf.calls)
    bf.offsets <- bf.onsets + sample(bdur, # and offset times
                                     size=bf.calls, replace=TRUE) 
    gf.offsets <- gf.onsets + sample(gdur, 
                                     size=gf.calls, replace=TRUE)
    species <- c(rep("B", bf.calls), rep("G", gf.calls))
    
    fd <- data.frame(species=species, onset=c(bf.onsets, gf.onsets), 
                     offset=c(bf.offsets, gf.offsets))
    fd$duration <- fd$offset - fd$onset
    
    fd <- getinterrupts(fd)
    return(p_gints(fd))
  }))
}

# load empirical data 
load("all_calls.RData")

ncalls <- commandArgs()[7]
ncalls <- as.integer(sub("\r", "", ncalls))
sims <- commandArgs()[8] 
print(commandArgs())
print(paste("ncalls = ", ncalls, "\nsims = ", sims, sep=""))
p.int <- fpow(data=all.calls, ncalls=ncalls, sims=sims)
pv.calc <- ecdf(p.int)
pv <- pv.calc(0)
int_sum <- list(ncalls=ncalls, nsims=sims, p.int=p.int, pv=pv)
save(int_sum, file=paste("./data/", ncalls, "calls.RData", sep=""))
