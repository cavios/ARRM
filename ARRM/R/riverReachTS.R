## Time-space model for altimetry data located on a river reach via an AR1 process
## in both the time and space directions


## Negative joint likelihood (nll) of data and parameters
f<-function(param){
    "[<-" <- ADoverload("[<-")
    getAll(data, param, warn = FALSE)
    phi<-2*plogis(logitphi)-1
    sd<-exp(logSd)
    sdEta<-exp(logSdEta)
    cumExpMu<-rep(0,length(mu))
    cumExpMu[1]=mu[1] 
    for(i in 2:length(mu))cumExpMu[i]=cumExpMu[i-1]+exp(mu[i])
    nll<-0
    sp<-splinefun(xknot,cumExpMu)
    spfield<-sp(dfield)
    v<-t(t(eta)-spfield)
    spObs<-sp(dist)
    deltasp<-spObs-spfield[xidx]
    

    df1 <- function(x)dautoreg(x, phi=phi[1],scale=sdEta,log=TRUE)
    df2 <- function(x)dautoreg(x, phi=phi[2], log=TRUE)
    nll<-nll - dseparable(df1, df2)(v)

    pred=eta[cbind(tidx,xidx)]+deltasp
    id<-which(satid < (length(bias)+1))
    if(length(id)!=0) pred[id]<-pred[id]+bias[satid[id]]
        
    predsd<-sd[satid]
    nll<-nll -sum(w*dnorm(H,pred,predsd,TRUE)) 

    residual<-(H-pred)/sd[satid]
    REPORT(pred)
    REPORT(predsd)
    REPORT(cumExpMu)
    REPORT(eta)
    
    
    tsVS=eta[,gid[1]]
    predsd<-sd[satid]
    ADREPORT(tsVS)
    nll
}

##' Fit riverReachTS
##' @param dat Data
##' @param VS Reach distance of virtual stations 
##' @param tdim Number of time steps in the model
##' @param xdim Number of space steps in the model
##' @param Nknot Number of spline nodes to describe the mean river reach profile
##' @param nit Number of model iterations, usually 2-5 is enough
##' @param cutp The fraction of expected outliers
##' @import RTMB
##' @examples
##' fit<-riverReachTS(dat)
##' @export

riverReachTS<-function(dat, VS=NULL,tdim=NULL,xdim=20,Nknot=10,nit=5,cutp=0.05 ){

#prepare data
    dis<- round(diff(range(dat$chainage)))
    dat$satid<-as.integer(as.factor(dat$satid))
    o<-order(dat$time)
    dato<-dat[o,]

    startH<-mean(dato$height[order(dato$chainage)[1:100]])
    if(is.null(tdim)){
        mtime<-round(mean(diff(unique(dato$time))*365),0)
        tdim<-round(diff(range(dat$time))*365/mtime,0)
    }
    if(is.null(VS)){
        VS<-diff(range(dato$dist))/2
    }
    

    tidx <- assignPretty(dato$time, tdim)
    xidx <- assignPretty(dato$chainage, xdim)
    gid <- assignPretty(VS, xdim)$ct
    x<-xidx$xx
    y<-tidx$xx


    myrange<-range(dato$chainage)
    xrange<-range(dato$chainage)/xdim
    xknot<-seq(myrange[1],myrange[2],length=Nknot)

    data<-list(H=dato$height,
           xidx=xidx$ct,
           tidx=tidx$ct,
           dfield=xidx$xx,
           satid=dato$satid,
           w=0*dato$height+1,
           dist=dato$chainage,
           xknot=xknot,
           VS=VS,
           gid=gid)



    param<-list(eta=matrix(0,tdim,xdim),
                logSd=rep(0,length(unique(data$satid))),
                logSdEta=0.1,
                logitphi=rep(invf(0.5),2),
                mu=c(startH,rep(0,length(data$xknot)-1)),
                bias=rep(0,length(unique(data$satid))-1)
                )
    
    data <- local({data<-data; environment()})
    environment(f) <- data
    
    
    allnll<-numeric(nit)
    obj <- MakeADFun(f,param,random=c("eta"), inner.control = list(maxit = 100),ridge.correct = FALSE,silent=TRUE)
    opt<-nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=5000, iter.max=5000))
    
    allnll[1]<-opt$objective
    environment()
    mypred<-obj$report()$pred
    mypredsd<-obj$report()$predsd
    for(i in 2:length(allnll)){
        cat('doing itteation ',i,'\n')
        tailp <- 1-pnorm(abs(dato$height-mypred), mean=0, sd=mypredsd)
        w <- ifelse(tailp<cutp, tailp/cutp, 1)
        eps <- 1.0e-16
        data$data$w= ifelse(w<eps, eps, w)
        environment(f) <- data
        obj <- MakeADFun(f,param,random="eta", inner.control = list(maxit = 100),silent=TRUE)
        opt<-nlminb(obj$par,obj$fn,obj$gr, control=list(eval.max=3000, iter.max=3000))
        allnll[i]<-opt$objective
        environment()
        mypred<-obj$report()$pred
        mypredsd<-obj$report()$predsd
        
    }

    rep<-sdreport(obj)
    environment()
    mypred<-obj$report()$pred
    time<-tidx$xx
    pl<-as.list(rep, "Est")
    plsd<-as.list(rep, "Std")
    
    ts<-as.list(rep,"Est",report=TRUE)$tsVS
    tssd<-as.list(rep,"Std",report=TRUE)$tsVS
    slope<-exp(pl$mu)
    slope[1]<-pl$mu[1]
    slope<-cumsum(slope)

    myts<-data.frame(time=y,H=ts,Hsd=tssd)
    out<-list(x=x,y=y,eta=pl$eta,etasd=plsd$eta,tsVS=myts)
    class(out)<-"ARRM"
    return(out)
}



##' Plot an object returned by the function riverReachTS()
##' @param fit Object returned by get.TS()
##' @param dat The raw water level data can be added to the plot. 
##' @param ... Additional argumants to plot
##' @importFrom fields image.plot
##' @import viridis
##' @keywords plot
##' @export

plot.ARRM<-function(fit,dat=NULL,doSave=FALSE, type="F",...){
    x<-fit$x
    y<-fit$y
    eta<-fit$eta
    tsvs<-fit$tsVS
    myrange<-range(c(tsvs$H+2*tsvs$Hsd,tsvs$H-2*tsvs$Hsd))
    if(doSave){pdf('myfield.pdf',12,10)}
    par(mfrow=c(2,1),mar=c(4,4,1,1))
    image.plot(y,x,eta,xlab='Time in decimal years',ylab='Reach centerline distance [m]',col=viridis::turbo(50))
    plot(tsvs$time,tsvs$H,t='l',col='blue',lwd=5,ylim=myrange,xlab='Time in decimal years',ylab='Elevations w.r.t. EGM2008')
    
    mypoly<-list()
    mypoly$x<-c(y,rev(y))
    mypoly$y<-c(tsvs$H-2*tsvs$Hsd,rev(tsvs$H+2*tsvs$Hsd))
    polygon(mypoly,col='azure2',border=NA)
    lines(tsvs$time,tsvs$H,col='blue4',lwd=3)
    #lines(tsvs$time,tsvs$H+2*tsvs$Hsd,t='l',col='gray',lwd=1,lty=1)
    #lines(tsvs$time,tsvs$H-2*tsvs$Hsd,t='l',col='gray',lwd=1,lty=1)
    if(!is.null(dat)){
        points(dat$time,dat$height,col='gray',pch=dat$satid)
        lines(tsvs$time,tsvs$H,col='blue4',lwd=5)
    }

    if(doSave){dev.off()}
}


##' Plot timeseries of an object returned by the function riverReachTS()
##' @param fit Object returned by get.TS()
##' @param dat The raw water level data can be added to the plot. 
##' @param ... Additional argumants to plot
##' @importFrom fields image.plot
##' @import viridis
##' @keywords plot
##' @export

TSplot<-function(fit,dat=NULL,...){
    x<-fit$x
    y<-fit$y
    eta<-fit$eta
    nc=length(x)
    matplot(y,eta,t='l',col=viridis::viridis(nc),lty=1,xlab='Time in decimal years',ylab='Surface elevation w.r.t. EGM2008',...)
    if(!is.null(dat)){
        points(dat$time,dat$height,col='gray',pch=dat$satid)
        
    }
}



plotProfile<-function(fit){
    x<-fit$x
    y<-fit$y
    eta<-fit$eta
    nc=length(x)
    matplot(x/1000,t(eta),t='l',col=viridis::viridis(nc),lty=1,xlab='chainage km',ylab='Surface elevation w.r.t. EGM2008')
}


plotDat<-function(dat){

    par(mfrow=c(2,2))
    mycol <- viridis::viridis(50)[as.numeric(cut(dat$height,breaks = 50))]
    plot(dat$lon,dat$lat,col=mycol,pch=16,xlab='Longitude',ylab='Latitude')
    plot(dat$chainage/1000,dat$height,col=dat$satid,xlab='Chainage km',ylab='Elevation w.r.t. EHM2008')
    plot(dat$time,dat$chainage/1000,col=mycol,xlab='Time in decimal years',ylab='Chainage km')
    
}


