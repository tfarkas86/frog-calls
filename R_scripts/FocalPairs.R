rm(list=ls())
library(parallel)

##### FUNCTIONS ###############################################################

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

rand <- function(td, bounded=TRUE, min=0, max=1800){
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

png(file=paste("./pdfs/", int_sum$id, ".png", sep=""))

den <- density(int_sum$rand.p.ints, from=0, adjust=2)
plot(x=den, xlab="Proportion Green Frog Calls Interrupting", type="l", 
     ylab="Probability Density", main=paste("id: ", int_sum$id, "; ", 
                                            "p = ", int_sum$pv, "; ", 
                                            "nboot = ", int_sum$enn, sep=""))
abline(v=int_sum$ci[1], lty=2)
arrows(0, max(den$y)/1.5, 0, 0, col="red", lwd=2)

dev.off()
