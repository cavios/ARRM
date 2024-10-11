

getProjUTM<-function(x,y){
    hem<-ifelse(mean(y)< 0,'south','north')
    myutm<-floor((mean(x) + 180) / 6) + 1
    myProj<-paste0("+proj=utm +zone=",myutm," +",hem," +ellps=WGS84")
    myProj
}

dist2<-function(v, w)(v[1] - w[1])^2 + (v[2] - w[2])^2
getdist<-function(v, w)sqrt((v[1] - w[1])^2 + (v[2] - w[2])^2)
getdist2<-function(v, w)sqrt((v[1] - w[,1])^2 + (v[2] - w[,2])^2)

distToSegmentSquared<-function(p, v, w) {
  l2 = dist2(v, w);
  if (l2 == 0){ 
    ret<-c(dist2(p, v),v);
  }else{
    t = ((p[1] - v[1]) * (w[1] - v[1]) + (p[2] - v[2]) * (w[2] - v[2])) / l2;
    if (t < 0) ret<- c(dist2(p, v),v);
    if (t > 1) ret<- c(dist2(p, w),w);
    if((t>=0)&(t<=1)){
      proj<-c(v[1] + t * (w[1] - v[1]),v[2] + t * (w[2] - v[2]) )
      ret<-c(dist2(p, proj),proj);
    }
  }
  return(ret)
}


doone2<-function(xydat,xy,mydistm1){
    if(mydistm1==0){
        eval<-sapply(1:1,function(i)distToSegmentSquared(xydat,xy[mydistm1+i,],xy[mydistm1+i+1,]))  
    }else if (mydistm1==nrow(xy)-1){
        eval<-sapply(0:0,function(i)distToSegmentSquared(xydat,xy[mydistm1+i,],xy[mydistm1+i+1,]))
    }else{ 
        eval<-sapply(0:1,function(i)distToSegmentSquared(xydat,xy[mydistm1+i,],xy[mydistm1+i+1,]))
    }
    
    pick<-which.min(eval[1,])
    c(eval[2:3,pick],pick)
}


##' Converts SWOT data to ARRD format.
##' @param swot RiverSP SWOT data
##' @param addChainage if TRUE chainage information will be added
##' @param CL if addChaimage=TRUE a centerline must be provided
##' @import lubridate 
##' @examples
##'out<-SWOT2AARD(swot,addChainage=TRUE,CL=CL)
##' @export

SWOT2AARD<-function(swot,addChainage=FALSE,CL=NULL){
    mydate<-substr(swot$time_str,1,10)
    date <- ymd(mydate)
    decYear<-decimal_date(date) # 2009.11
    dat<-cbind(swot,decYear)

    if(addChainage)dat<-addChainage(dat,CL)
    o<-order(dat$decYear)
    dat<-dat[o,]
    out<-data.frame(timesec=dat$time, time=dat$decYear,  lat=dat$lat, lon=dat$lon, height=dat$wse, distnode=NA, width=dat$p_width, reachID=dat$reach_id, nodeID=dat$node_id, wse=dat$p_wse, OCval=NA, satid=12)
    if(addChainage){
        chainage<-dat$chainage
        out<-cbind(out,chainage)
    }
    out
}

##' Project altimetry data to the river centerline.
##' access to a basin sqlite file. Basin sqlite files can be acessed from ...  
##' @param dat an oject of the function getAltReaches or a data.frame that contains at least a column with latitude and longitude named "lat" and "lon", respectively.  
##' @param centerline an object from the function extractCL or matrix/dataframe with centerline coordinates of latitude and longitude named "y" and "x", respectively.     
##' @importFrom terra vect
##' @importFrom terra project
##' @importFrom terra crds
##' @examples
##'out<-projectToCenterline(dat,centerline)
##' @export

addChainage<-function(dat,centerline){
    datLL<-cbind(dat$lon,dat$lat)
    centerLineLL<-cbind(centerline$x,centerline$y)
   #transform data to UTM coordinates
    myProj<-getProjUTM(centerLineLL[,1], centerLineLL[,2])
    v <- vect(datLL, crs="+proj=longlat")
    u <- terra::project(v, myProj)
    xydat <- crds(u)

    vv <- vect(as.matrix(centerLineLL), crs="+proj=longlat")
    uu <- terra::project(vv, myProj)
    xy2 <- crds(uu)
    xy<-xy2[!duplicated(xy2),]
    
    
# calculate length of line segments of centerLine
    NP<-nrow(xy)
    BeginPx<-xy[1:(NP-1),1]
    BeginPy<-xy[1:(NP-1),2]
    begin<-cbind(BeginPx,BeginPy)
    EndPx<-xy[2:NP,1]
    EndPy<-xy[2:NP,2]
    end<-cbind(EndPx,EndPy)

    lineSegLength<-sqrt((EndPx-BeginPx)^2+(EndPy-BeginPy)^2) 
    alpha<-(EndPy-BeginPy)/(EndPx-BeginPx)
#project points to lines
    
    mydist<-unlist(lapply(1:nrow(xydat),function(i)which.min(getdist2(xydat[i,],xy))))
    mydistm1<-mydist-1
    all<-lapply(1:nrow(xydat),function(i)doone2(xydat[i,],xy,mydistm1[i]))
    all<-do.call(rbind,all)
    all[,3]<-all[,3]-1+mydistm1

    #identifi line segment and get dist
    MyPos<-matrix(,nrow=nrow(xydat),ncol=5)
    MyPos[,1:2]<-all[,1:2]
    yyy <- vect(MyPos[,1:2,drop=FALSE], crs=myProj)
    myll <- terra::project(yyy, "+proj=longlat +datum=WGS84")
    ll<-crds(myll)

    MyPos[,3:4]<-ll

    
    #smalldist<-rep(NA,nrow(all)) 
    for (i in 1:nrow(MyPos)){
        
        id<-which.min(sqrt((all[i,1]-begin[,1])^2+(all[i,2]-begin[,2])^2))

        
        if (abs(all[i,1]-begin[id,1])<0.00001){
            alphap<-(all[i,2]-end[id,2])/(all[i,1]-end[id,1])
        } else{alphap<-(all[i,2]-begin[id,2])/(all[i,1]-begin[id,1])}
        
        myna<-sum(is.na(alpha))
        test<-abs(alphap-alpha[id])<0.0000000001
        
        if (test){
            if (id==1) distTOT<-getdist(all[i,],begin[id,])
            else distTOT<-sum(lineSegLength[1:(id-1)])+getdist(all[i,],begin[id,])
        } else {
            if (id==2) distTOT<-getdist(all[i,],begin[id-1,])
            else if(id==1) distTOT<-getdist(all[i,],begin[id,])
            else distTOT<-sum(lineSegLength[1:(id-1-1)])+getdist(all[i,],begin[id-1,])
        }
        MyPos[i,5]<-distTOT
        
    }
    chainage<-MyPos[,5]
    out<-cbind(dat,chainage)
    return(out)
}





getOneReach<-function(mydb,myreach,mybasin){
    myq<-paste0("SELECT * FROM basin",mybasin," WHERE reachID =" ,myreach)
    res <- dbSendQuery(mydb, myq)
    out<-dbFetch(res)
    dbClearResult(res)
    out
}

##' Extract nadir altimetry data from the Altimetry River Reach database (ARRD). This function requires
##' access to a basin sqlite file. Basin sqlite files can be acessed from ...  
##' @param dbfile a particular basin sqlite file e.g. "7428.sqlite". The basin number is the HydroShed level 4 PFAF_ID 
##' @param mybasin HydroShed basin PFAF_ID (level 4) 
##' @param myReaches a vector with SWORD reach numbers
##' @param addFilter =FALSE, if TRUE only observation with a node distanstance smaller than the width of the river is returned
##' @param addBuffer =NULL, if a value is specified only observation with a node distanstance smaller than the width of the river + (-) addBuffer is returned
##' @param removeNA = TRUE, remove rows where the elevation is NA
##' @import DBI
##' @import RSQLite
##' @examples
##' dbfile<-"7428.sqlite"
##' mybasin<-'7428'
##' myreaches<-c(74287400183, 74287400163, 74287400151, 74287400173)
##' out<-getAltReaches(dbfile,mybasin,myreaches)
##' @export

getAltReaches<-function(dbfile,mybasin,myreaches,addFilter=FALSE,addBuffer=NULL,removeNA=TRUE){
    mydb <- dbConnect(RSQLite::SQLite(), dbfile)
    out<-lapply(1:length(myreaches),function(i){cat('doing ',i,'\n');getOneReach(mydb,myreaches[i],mybasin)})
    out<-do.call(rbind,out)
    dbDisconnect(mydb)
    if(removeNA){
        id<-which(!is.na(out$height))
        out<-out[id,]
    }
    if(addFilter){
        id<-which(out$distnode < out$width)
        if(!is.null(addBuffer)){
            id<-which(out$distnode < out$width+addBuffer)
        }
        out<-out[id,]
    }
    out
}










