
setwd("C:/Users/Susan Herrick/Research/Acoustics Chapter/R work for Acoustics/Analysis")

##############################################
### now skip to plotting section below

library(parallel);library(abind);library(reshape)

### read in files
files=list.files(path="~/Dropbox/Projects/FrogCalls/data", 
                 pattern="1.csv" , full.names=TRUE)  
#get filenames

getdata<-function(i) {
  id=substr(i,6,13) 
  dtemp=read.csv(i) 
  gf=dtemp[,c("file","onset","duration")]
  bf=dtemp[,c("file.1","onset.1","duration.1")];colnames(bf)=c("file","onset","duration")
  dtemp=rbind(gf,bf)
  dtemp=dtemp[!is.na(dtemp$duration),]
  print(i)
  return(dtemp)
}


### Import data
d=do.call(rbind.data.frame,lapply(files,getdata))

d <- fd
### Organize and add a few variables
d$id=factor(substr(d$file,2,8))
d$species=factor(substr(d$file,1,1),label=c("Green Frog","Bullfrog"))
d$month=as.numeric(substr(d$file,3,3))
d$day=as.numeric(substr(d$file,4,5))

d$minute=as.numeric(substr(d$onset,1,2))
d$second=as.numeric(substr(d$onset,4,7))

## Drop weird records with AM in time
d=d[!is.na(d$second),]

## calculate onset and offset
d$onset=d$minute*60+d$second
d$offset=d$onset+d$duration

d=d[d$onset<=1800,]

### Add flag if call starts during call of other species
  getinterrupts=function(td,random=F){
  td2=do.call(rbind.data.frame,lapply(1:nrow(td),
    function(l){
    sp=td$species[l]
    tonset=td$onset[l]
    w=tonset>=td$onset&tonset<=td$offset  #find overlapping calls  (including self)
    w[l]=F  #set self to false
    w=which(w)[1]  #get index of any trues
    ## get most recent but ended offset time
    b=which(td$offset<=tonset)  #find all previously ending calls
    spid=td$species[b]==td$species[l]
    ptimehom=(tonset-td$offset[b][spid])[which.min(tonset-td$offset[b][spid])]
    ptimehet=(tonset-td$offset[b][!spid])[which.min(tonset-td$offset[b][!spid])]
    if(length(w)>1) print("interrupted more than one call")
    ttd=data.frame(td[l,],  #drop these cols if present
      interrupt=ifelse(length(w)==0,"none",as.character(td$species[w])),
      itime=ifelse(length(w)==0,NA,tonset-td$onset[w]),
      otime=ifelse(length(w)==0,NA,min(td$offset[w],td$offset[l])-tonset),
      ptimehet=ifelse(length(ptimehet)>0,ptimehet,NA),ptimehom=ifelse(length(ptimehom)>0,ptimehom,NA))
    return(ttd)
  }))
print(as.character(td$id[1]))
return(td2)}

### Run function on real data
d=do.call(rbind.data.frame,lapply(split(d,d$id),getinterrupts,random=F))

### add species code
d$speciescode=ifelse(d$species=="Bullfrog","bf","gf")

### add minute bin for plotting
d$bin=round(d$onset/60)

## create long format for plotting
dl=melt.data.frame(d,id.vars=c("id","species","onsetmin"),measure.vars=c("duration","itime","otime","ptimehet","ptimehom"))
dl=dl[!is.na(dl$value),]

durationsum=tapply(d$duration,list(d$id,d$speciescode),sum,na.rm=T)

### Summarize variables

#    otime=100*tapply(x$otime,split,sum,na.rm=T)/1800,   #old way
#    otime=100*tapply(x$otime,tsplit,sum,na.rm=T)/tapply(x$duration,split,sum,na.rm=T),

summaryfun=function(x,split=c("id","speciescode"),bin=F,durationtable=durationsum) {
  tsplit=as.list(x[,split])
  otime=100*tapply(x$otime,tsplit,sum,na.rm=T)/durationsum
  summary=data.frame(
  otime=otime,
     itime=tapply(x$itime,tsplit,mean,na.rm=T),
    ptimehet=tapply(x$ptimehet,tsplit,mean,na.rm=T),
    ptimehom=tapply(x$ptimehom,tsplit,mean,na.rm=T),
    gfbfinterrupt=100*tapply(ifelse(x$species=="Green Frog"&x$interrupt=="Bullfrog",1,0),tsplit,sum,na.rm=T)/sum(x$species=="Green Frog"),
    gfgfinterrupt=100*tapply(ifelse(x$species=="Green Frog"&x$interrupt=="Green Frog",1,0),tsplit,sum,na.rm=T)/sum(x$species=="Green Frog"))
  if(any(grepl("duration",colnames(x)))) summary[,c("pduration.bf","pduration.gf")]=100*tapply(x$duration,tsplit,sum,na.rm=T)/1800
  summary$id=rownames(summary)
  summary2=melt(summary,id.var="id")
  t=do.call(rbind,strsplit(as.character(summary2$variable),"[.]"))
  summary2$species=t[,2]
  if(bin) summary2$bin=t[,3]
  summary2$variable=t[,1]
  return(summary2)
}

ds=summaryfun(d,bin=F)
dsc=cast(ds,id+species~variable,value="value")

### summarize by minute
durationsummin=tapply(d$duration,list(d$id,d$speciescode,d$bin),sum,na.rm=T) ##############
dsm=summaryfun(d,split=c("id","speciescode","bin"),bin=T,durationtable=durationsummin)

###############################################################
### randomizations of onset times
n=2500
uid=unique(d$id)

randomize=function(td){
  tgf=td$species=="Green Frog"
  td$onset[tgf]=runif(sum(tgf),0,1800)  ## shuffle the onset time
  td$offset[tgf]=td$onset[tgf]+td$duration[tgf]  #calculate new offset
  tr=getinterrupts(td)  #get indicies
  tr=tr[,c("file","id","species","interrupt","itime","otime","ptimehet","ptimehom")]  #keep only a few columns
  return(tr)
  }

### run it and save
system.time(r<<-lapply(1:n,function(q) do.call(rbind.data.frame,mclapply(split(d[,-grep("interrupt|itime|otime|ptimehet|ptimehom",colnames(d))],d$id),
                                                                         randomize,mc.cores=4,mc.preschedule=F))))
save(r,file="output/randomizations.Rdata")

### Run the summary function on the randomizations 
rs=abind(lapply(r,summaryfun,split=c("id","speciescode"),bin=F),along=3)

### Melt summary to 'long' format
rsl=melt.array(rs,varnames=c("row","element","iteration"))
rsl=cast(rsl,row+iteration~element,value="value")
rsl$value=as.numeric(as.character(rsl$value))
rsl=rsl[!(rsl$species=="bf"&rsl$variable%in%c("gfbfinterrupt","gfgfinterrupt")),]
rsl$speciescode=rsl$species
rsl$species=as.factor(ifelse(rsl$species=="bf","Bullfrog","Green Frog"))

rst=as.data.frame(t(sapply(split(rsl$value,list(rsl$id,rsl$variable,rsl$species)),quantile,c(0.025,.5,.975),na.rm=T)))
rst[,c("id","variable","species")]=do.call(rbind.data.frame,strsplit(rownames(rst),".",fixed=T))
rownames(rst)=1:nrow(rst)

## table of interruptions
table(d$species,d$interrupt)

### output tables
write.csv(d,file="output/data.csv",row.names=F)
write.csv(dl,file="output/datalong.csv",row.names=F)
write.csv(ds,file="output/SummaryObserved1.csv",row.names=F)
write.csv(dsc,file="output/SummaryObserved2.csv",row.names=F)
write.csv(dsm,file="output/SummaryObservedMinuteByMinute.csv",row.names=F)
write.csv(rsl,file="output/randomizations.csv",row.names=F)
write.csv(rst,file="output/SummaryRandomizations.csv",row.names=F)


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
#### some plots

### read in data
d=read.csv(file="output/data.csv")
dl=read.csv(file="output/datalong.csv")
ds=read.csv(file="output/SummaryObserved1.csv")
dsc=read.csv(file="output/SummaryObserved2.csv")
dsm=read.csv(file="output/SummaryObservedMinuteByMinute.csv")
rsl=read.csv(file="output/randomizations.csv")
rsl$id=as.factor(rsl$id)
rst=read.csv(file="output/SummaryRandomizations.csv")

### load libraries
library(lattice);library(latticeExtra)

nid=length(unique(rst$id))  # get number of groups for use later

### write a pdf
pdf("Summary.pdf",width=11,height=8.5,useDingbats=F)
trellis.par.set(strip.background =list(col=grey(.9)),fontsize =list(text=14))  #set colors for lattice plots


## Scatterplot of gf vs. bf % calling times
dsmd=cast(dsm[dsm$variable=="pduration",],id+bin~species)
xyplot(gf~bf,groups=dsmd$id,data=dsmd,ylab="Green Frog % time calling",xlab="Bullfrog % time calling",sub="Minute by minute percentages",col=1:nid,pch=1:nid,key=list(space="right",text=list(as.character(unique(dsmd$id))),points=list(pch=1:nid,col=1:nid)),
       main="Frog Activity")

## calling times
boxplot(duration~species,data=fd,log="y",las=1,col="grey",ylab="Call Duration in seconds (log scale)",cex=.5,main="Comparison of Calling Duration")

## Metric distributions
p1=densityplot(~value|variable,data=dl,pch="",auto.key=T,panel=function(x,subscripts){
  tdl=dl[subscripts,]
  panel.densityplot(tdl$value[tdl$species=="Green Frog"],subscripts=subscripts,type=NA,col="green",darg=list(n=1000))
  panel.densityplot(tdl$value[tdl$species=="Bullfrog"],subscripts=subscripts,type=NA,col="blue",darg=list(n=50))
},main=paste("Distribution of observed values for each species"),scale=list(y=list(draw=F),relation="free",alternating=1),
                  xlab="Seconds",layout=c(1,5))
update(p1,ylim=lapply(p1$y.limits,function(x) c(x[1],x[2]*1.5)))


### Randomizations - density plots with arrows
for(v in unique(rsl$variable)){
  trange=pretty(c(rsl$value[rsl$variable==v],ds$value[ds$variable==v]))
  print(useOuterStrips(densityplot(~value|id+species,data=rsl[rsl$variable==v,],pch="",auto.key=T,panel=function(x,subscripts){
    trb=rsl[rsl$variable==v,][subscripts,]
    tid=trb$id[1]
    tvar=trb$variable[1]
    tsp=trb$speciescode[1]
    tobs=ds$value[ds$id==tid&ds$species==tsp&ds$variable%in%tvar]
    tquant=quantile(na.omit(x),c(0.025,0.975))
    panel.arrows(tobs,1000,tobs,0,col=ifelse(tobs<tquant[1]|tobs>tquant[2],"red","blue"),length=.1)
    panel.densityplot(na.omit(x),subscripts=subscripts,type=NA,col="black")
  },type="count",main=paste(v," (",max(rsl$iteration)," randomizations)"),scale=list(y=list(draw=F),alternating=1),xlim=c(min(trange)-.5,trange),asp=1,
                                   xlab="Arrow indicates observed value, density indicates random expectation",
                                   sub=expression(paste("A red arrow indicates observation is significantly different from 0 (  ",alpha==0.05,")")))))
}

## timeseries as dashes
mbin=180
onsetbin=cut(d$onset,seq(0,30*60,mbin),labels=paste(seq(0,29*60,mbin)/60,"-",seq(mbin,30*60,mbin)/60," minutes"))

xyplot(species~onset|onsetbin+id,data=d,type="l",panel=function(x,y,subscripts){
    i=unique(d[subscripts,"id"])
    t=unique(onsetbin[subscripts])
    dbf=d[d$species=="Bull Frog"&d$id==i&onsetbin==t,]
    panel.rect(dbf$onset,.75,dbf$offset,1.25,lwd=2,col="dodgerblue",border=NA)
    dgf=d[d$species=="Green Frog"&d$id==i&onsetbin==t,]
    panel.rect(dgf$onset,1.75,dgf$offset,2.25,lwd=2,col="green",border=NA)
    if(which.packet()[1]==1) panel.text(0,2.2,i,cex=.8,col="grey")
  },subscripts=T,strip=F,scales=list(x=list(relation="sliced",draw=F),alternating=1),layout=c(1,10,10),as.table=T,ylab="Species",xlab="Time (3 minutes)"
)

### make colored timeseries
#levelplot(frog~timemin*id,data=o,col.regions=c("white","red","green","black"),ylab="Observation",xlab="Time (min)",
#          colorkey=list(at=c(0:4),labels=list(at=c(.5,1.5,2.5,3.5),labels=as.character(unique(o$frogtype))),height=.3),
#          main="Timeseries of calls",asp=.2)


dev.off()  #stop writing to the pdf
system("ps2pdf Summary.pdf Summary_small.pdf")  #shrink the file (only works if you have ps2pdf installed

