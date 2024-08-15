


assignPretty <- function(x,n){
  pr <- seq(min(x)-1/365,max(x)+1/365,length=n+1)
  mydiff<-diff(pr)[1]/2
  etatime<-pr[1:n]+mydiff
  br <- pr
  ct <- as.integer(cut(x,br))
  return(list(ct=ct, xx=etatime, br=br))
}


invf <- function(y) -0.5 * log(2/(y + 1) - 1)
