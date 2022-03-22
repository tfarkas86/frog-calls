###############################################################################
# R script for analysis of nearest-neighbor frog calls, Herrick et al. 2015
# Author: TE Farkas, timothy.farkas@gmail.com
###############################################################################

# SECTION 1: For use with Sun Grid Engine. Reads raw data files (CSV), preps
# them for analysis, and saves bootstrap objects to files (n = number of CSVs)
# as .RData for post processing (Section 2). Run as array job with 10000 
# simulations per CSV (n = 24)

##### FUNCTIONS ######

# getdata() takes a path to a data file and imports and preps it for analysis

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

# rand() function takes a data set, parses up the total recording duration
# into units of maximum green frog calling duration, randomly assigns green 
# frog onset times, recalculates offset times using empirical calling  duration, 
# and checks for interrupts using getinterrupts()

rand <- function(td, bounded=TRUE, min=0, max=1800) {
  start <- ifelse(bounded, min(td$onset), min)
  end <- ifelse(bounded, max(td$offset), max)
  gf <- td$species=="G" # is green frog?
  max.dur <- max(td$duration[td$species=="G"])
  onset.pts <- seq(start, end, by=.3)
  td$onset[gf] <- sample(x=onset.pts, # draw gf onset time onset.pts
                         size=sum(gf), replace=FALSE) 
  td$offset[gf] <- td$onset[gf]+td$duration[gf] # calculate new offset
  #tr <- getinterrupts(td) # get interrupts
  #tr <- tr[,c("species", "int_sp", "onset")] # keep only a few columns
  #return(tr)
  return(td)
}

# boot() function takes a data set, passes it to rand() and iterates enn times
# returning a list of randomized data sets

boot <- function(fd, enn, ...) {
  rands <- mclapply(1:enn, function(q) rand(fd, ...), mc.cores=8, 
                    mc.preschedule=FALSE)
  return(rands)
}

# getinterrupts() function takes a file, modified by getdata(), and adds 
# information as to which calls interrupt other calls, and the species 
# information for those interruptions

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

# p_gints() function gets proportion green frog calls that interrupt a bullfrog
# call and returns the proportion (p_ints)

p_gints <- function(td) {
  n_gcalls <- sum(td$species=="G")
  n_gints <- sum(td$species=="G" & td$int_sp=="B")
  p_gints <- n_gints/n_gcalls
  return(p_gints)
}

# n_bcalls() function returns number of bullfrog calls 

n_bcalls <- function(td) {
  n_bcalls <- sum(td$species=="B")
  return(n_bcalls)
}

# n_gcalls() function return number of green frog calls

n_gcalls <- function(td) {
  n_gcalls <- sum(td$species=="G")
  return(n_gcalls)
}

# n_gints function returns number of green frog interruptions

n_gints <- function(td) {
  n_gints <- sum(td$species=="G" & td$int_sp=="B")
  return(n_gints)
}

# ht.stats (hypothesis testing stats) is a function that calculates confidence 
# intervals and p-values for a vector based on an empirical cumulative density 
# function, and calculates depending on what alpha is desired, and whether 
# p-values are to be one or two sided, upper or lower

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

# int.stats() is a function calculating interruption statistics that
# returns a list with frequency of interrupts for randomized data ($rand_p_ints)
# and frequency for empirical data ($emp_p_ints). Pass arguments for rand() 
# including bounded(T/F) and min/max, and boot (enn = number of bootstraps)

int.stats <- function(file, enn=10, ...) {
  fd <- getdata(file) # reshape data
  id <- fd$id[1]
  rd <- boot(fd=fd, enn=enn, ...) # get list of randomized data sets
  r.ints <- lapply(rd, getinterrupts) # add interruption data to each rand
  e.ints <- getinterrupts(fd) # add interruption data to data file
  n.bcalls <- n_bcalls(fd) # get number bull calls
  n.gcalls <- n_gcalls(fd) # get number of green frog calls
  n.r.gints <- sapply(r.ints, n_gints) # get number of gf ints for rands
  p.r.ints<- sapply(r.ints, p_gints) # get proportion interrupts for rands
  n.e.gints <- n_gints(e.ints) # get number of gf ints for data
  p.e.ints <- p_gints(e.ints) # get proportion interrupts for data
  call.stats <- list(id=as.character(id), n.bcalls=n.bcalls, n.gcalls=n.gcalls, 
                     emp.n.ints=n.e.gints, emp.p.ints=p.e.ints, enn=enn, 
                     rand.n.ints=n.r.gints, rand.p.ints=p.r.ints)
  ht.stats <- ht.stats(data=p.r.ints, # get confidence and p-values
                       emp=p.e.ints, ...) 
  return(append(call.stats, ht.stats))
}

##### IMPORT AND ANALYSIS #####

file <- paste("./data/", commandArgs()[7], sep="")

int_sum <- int.stats(file=file, enn=10)
save(int_sum, file=paste("./stats/", int_sum$id, ".RData", sep=""))

# Section 2: Iteratively loads .RData files and analyzes int_sum objects for 
# each raw data file.

# get names of all int_sum data files
files <- list.files("~/Dropbox/Projects/FrogCalls/data/stats 2", full=TRUE)

# get mean % interruptings for each recording
avg.p.ints <- sapply(1:length(files), function(x) {
  load(files[x])
  return(mean(int_sum$rand.p.ints))
})

# get average number of interruptings (rounded down) 
avg.n.ints <- sapply(1:length(files), function(x) {
  load(files[x])
  return(floor(mean(int_sum$rand.n.ints)))
})

# get total number of green frog calls 
n.gcalls <- sapply(1:length(files), function(x) {
  load(files[x])
  return(int_sum$n.gcalls)
})

# get total number of bullfrog calls
n.bcalls <- sapply(1:length(files), function(x) {
  load(files[x])
  return(int_sum$n.bcalls)
})

ncalls <- n.gcalls + n.bcalls

##### MERGED DENSITY FUNCTION #####

# get % calls interrupting for all simulations for all raw data files
thsh <- 237 # minimum number of total calls (see Section 3: Power Analysis)
rands <- do.call(rbind, lapply(1:length(files), function(x) {
  load(files[x])
  if(int_sum$n.bcalls + int_sum$n.gcalls > 237) {
    return(100*int_sum$rand.p.int)
  }
}))

# plot merged histogram
par(mar=c(5,5,2,2))
conf <- quantile(rands, c(.05))
hist(rands, breaks=30, xlim=c(0, 30), xlab="Percent of calls interrupting", 
     main="", yaxt="n", asp=1, ylim=c(0, 25000), cex.axis=.75, cex.lab=.75,
     ylab="Number of simulations / 1000")
axis(2, at=seq(0, 25000, by=5000), labels=seq(0, 25, by=5), las=1, cex.axis=.75)
lines(x=rep(conf,2), y=c(0, 21000), lty=2, lwd=2)
text(x=conf-1, y=19500, labels="95% CI", cex=.75, srt=90)

pvs <- ecdf(rands) # create merged cumulative empirical density function
pvs(0) # get p-value (proportion 0s from )

# Section 3: Power analysis to determine minimum level of calling activity at 
# which significant interruption avoidance can be detected at alpha = 0.05.

# get all raw data files
emp.files <- list.files("~/Dropbox/Projects/FrogCalls/data", pattern=".csv",
                        full=TRUE)

# make huge dataframe of all the calls
all.calls <- do.call(rbind, lapply(emp.files, function(x) {
  getdata(x) 
}))
all.calls <- all.calls[-which(all.calls$id=="0624051"),]
all.calls$id <- as.factor(as.character(all.calls$id))
all.calls <- all.calls[all.calls$duration < 2,] # remove bogus bf call 
save(all.calls, file="~/Dropbox/Projects/FrogCalls/power/all_calls.RData")

# fpow() takes a data frame with all calling data, the total number of calls 
# in a recording session, and the number of simulations to be performed. It 
# returns the proportion of green frog calls interrupting bullfrog calls. For 
# use with Sun Grid Engine as an array job across calling activity (50:600, 10)

fpow <- function(all.calls, ncalls, sims=10) {
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

ncalls <- commandArgs()[7] # take ncalls as argument from bash script
sims <- commandArgs()[8] # take sims as argument from bash script

p.int <- fpow(data=all.calls, ncalls=100, sims=20)
pv.calc <- ecdf(p.int) # make emprirical cumulative density function
pv <- pv.calc(0) # calculate p-value
int_sum <- list(ncalls=ncalls, nsims=sims, p.int=p.int, pv=pv) # make object
save(int_sum, file=paste("./data/", ncalls, "calls.RData", sep="")) # save

# get list of files from power analysis
files <- list.files("~/Dropbox/Projects/FrogCalls/power/data.10000/", 
                    full=TRUE)

# get vector of calling number
ncalls <- do.call(c,(lapply(files, function(x) {
  load(x)
  int_sum$ncalls
})))

# get vector of p values
pvs <- do.call(c,(lapply(files, function(x) {
  load(x)
  int_sum$pv
})))

# make data from with ncalls and pvs
cd <- data.frame(ncalls, pvs, sims)
cd <- cd[order(cd$ncalls),] # sort by ncalls

approx(x=cd$pvs[19:20], y=cd$ncalls[19:20], xout=.05) # interpolate = 237
plot(cd$ncalls, cd$pvs, type="l")
abline(h=.05, lty=2)

