library(foreach)
library(doParallel)

test = c()

for(i in 1:10){
  x = rnorm(1000)
  y = rnorm(1000)
  
  test[i] = lm(x~y)$coef[2]
}

  cl = makeCluster(1)
  registerDoParallel(cl)

  test.par = foreach(i = 1:10)%dopar%{
    
  
      
                      x = rnorm(1000)
                      y = rnorm(1000)
                      
                      test[i] = lm(x~y)$coef
                    }
  stopCluster(cl)