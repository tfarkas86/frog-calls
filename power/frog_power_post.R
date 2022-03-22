rm(list=ls())

files <- list.files("~/Dropbox/Projects/FrogCalls/power/data.10000/", 
                    full=TRUE)

load(files[1])

ncalls <- do.call(c,(lapply(files, function(x) {
  load(x)
  int_sum$ncalls
})))

pvs <- do.call(c,(lapply(files, function(x) {
  load(x)
  int_sum$pv
})))

sims <- do.call(c,(lapply(files, function(x) {
  load(x)
  int_sum$nsims
})))

cd <- data.frame(ncalls, pvs, sims)
cd <- cd[order(cd$ncalls),]

an1 <- lm(ncalls ~ exp(pvs), data=cd)
newdata <- data.frame(pvs=.05)
predict.lm(an1, newdata=newdata)
approx(x=cd$pvs[19:20], y=cd$ncalls[19:20], xout=.05)
plot(cd$ncalls, cd$pvs, type="l")
abline(h=.05, lty=2)
