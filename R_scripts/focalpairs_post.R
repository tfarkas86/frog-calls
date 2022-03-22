rm(list=ls())

# get names of all int_sum data files
files <- list.files("~/Dropbox/Projects/FrogCalls/data/stats 2", full=TRUE)

# get mean % interruptings for each recording
avg.p.ints <- sapply(1:length(files), function(x) {
  load(files[x])
  return(mean(int_sum$rand.p.ints))
})

##### one-sample t-test

ang.avgs <- asin(sqrt(avg.p.ints)) # angular transformation

tee <- t.test(x=ang.avgs, alternative="two.sided") # one-sample t-test (vs. 0)

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

adj.pees <- p.adjust(pees, method="holm")
min(ncalls[which(adj.pees < 0.05)])
adj.pees[which(ncalls > 250)]


avg.n.xint <- n.gcalls - avg.n.ints # number of gf calls NOT interrupting
int.dep <- cbind(avg.n.ints, avg.n.xint) # make 2 column dependent variable
s.int.dep <- int.dep
s.int.dep[,1] <- s.int.dep[,1] + floor(s.int.dep[,2]/2) # round down = conserve
s.int.dep[,2] <- s.int.dep[,2] - ceiling(s.int.dep[,2]/2)

### BINOMIAL GLM

an1 <- glm(s.int.dep ~ 1, family=quasibinomial) # OD = 1.08
summary(an1) # b = -2.75 (0.060 prob), se = 0.087, z = 32.94

##### FISHERS COMBINED P-VALUE METHOD ######

# get p-values
pees <- sapply(1:length(files), function(x) {
  load(files[x])
  return(int_sum$pv)
})

Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- 1-pchisq(Xsq, df = 2*length(p))
  return(c(Xsq = Xsq, p.value = p.val))
}

library(MADAM)
pees <- matrix(pees, nrow=1)
pees <- pees + .048 # bonferroni correction
fisher.method(pees) # p = 7.7 x 10e-8

## MERGED DENSITY FUNCTION ##

rands <- do.call(rbind, lapply(1:length(files), function(x) {
  load(files[x])
  if(int_sum$n.bcalls + int_sum$n.gcalls > 237) {
  return(100*int_sum$rand.p.int)
  }
}))

par(mar=c(5,5,2,2))
conf <- quantile(rands, c(.05))
bks <- -1:33
h1 <- hist(rands, breaks=bks, xlim=c(0, 30), xlab="Percent of calls interrupting", 
     main="", yaxt="n", asp=1, ylim=c(0, 25000), cex.axis=1, cex.lab=1,
     ylab=expression(paste("Number of simulations * 10"^"-3")))
axis(2, at=seq(0, 25000, by=5000), labels=seq(0, 25, by=5), las=1, cex.axis=1)
lines(x=rep(conf,2), y=c(0, 21000), lty=2, lwd=2)
text(x=conf-1, y=19500, labels="95%", cex=1, srt=90)
arrows(x1=h1$mids[1], y1=h1$counts[1], x0=h1$mids[1], y0=10000, angle=30,
       length=.15, lwd=2)

plot(density(rands, adjust=2.5, from=0))
pvs <- ecdf(rands)

pvs(0) # 0.004

rands <- as.vector(rands)
library(lattice)
densityplot(x=~rands, type="count", asp=1, 
            scale=list(y=list(draw=F), alternating=1), n=10)



## HISTOGRAMS ###

pdf("~/Dropbox/Projects/FrogCalls/pdfs/composite.pdf")

par(mfrow=c(6, 4), mar=c(2.5, 2, 2, 1)) # make 4 x 6 composite panel

lapply(1:length(files), function(i) {
  load(files[i])
  ncalls <- int_sum$n.bcalls + int_sum$n.gcalls
  if (ncalls > 281) {
  pv <- p.adjust(int_sum$pv, method="holm", n=1)
  l.ci <- quantile(int_sum$rand.p.ints, probs=0.05/1)
  hist(int_sum$rand.p.ints, 
       main=NULL, ylab=NULL, xlab=NULL)
  title(paste("id: ", int_sum$id, "; p = ", pv, ";", "\n", "ncalls = ", 
              int_sum$n.bcalls + int_sum$n.gcalls, sep=""), cex.main=1)
  abline(v=l.ci, lwd=2, lty=2) 
  }
})

dev.off()

##### calculate ncalls cutoff #####

# functions
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
#

emp.files <- list.files("~/Dropbox/Projects/FrogCalls/data", pattern=".csv",
                        full=TRUE)
emp.files <- emp.files[-length(emp.files)]

# make huge list of all the calls
all.calls <- do.call(rbind, lapply(emp.files, function(x) {
  getdata(x)
})
)

all.calls <- all.calls[-1183,] # remove bogus bf call (20 seconds!)

bdur <- all.calls[all.calls$species=="B", "duration"] # get all bf durations
gdur <- all.calls[all.calls$species=="G", "duration"] # get all gf durations

gprop <- n.gcalls/ncalls # proportion of total calls by green frogs
onset.pts <- seq(0, 1798.5, by=.3) # get set of onset times

# function computes frequency distrubutions for interruptions at variety
# of calling activity levels (from 50 to 600 calls in an hour). Next, consider 
# increasing number of simulations from 10 to 1000, and making an array job 
# on SGE

mclapply(seq(50, 600, 10), mc.cores=4, function(enn.calls) {
  do.call(c, mclapply(1:10, function(x) {
    gf.calls <- round(sample(g.prop, size=1)*enn.calls) # number of gf calls
    bf.calls <- enn.calls - gf.calls # number of bf calls
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
})

