## simple script joins all frog call data from 24 files

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

emp.files <- list.files("~/Dropbox/Projects/FrogCalls/data", pattern=".csv",
                        full=TRUE)
# make huge data.frame of all the calls
all.calls <- do.call(rbind, lapply(emp.files, function(x) {
  getdata(x) 
}))
all.calls <- all.calls[-which(all.calls$id=="0624051"),]
all.calls$id <- as.factor(as.character(all.calls$id))
all.calls <- all.calls[all.calls$duration < 2,] # remove bogus bf call 
save(all.calls, file="~/Dropbox/Projects/FrogCalls/power/all_calls.RData")