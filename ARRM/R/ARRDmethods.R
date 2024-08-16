

getProjUTM<-function(x,y){
    hem<-ifelse(mean(y)< 0,'south','north')
    myutm<-floor((mean(x) + 180) / 6) + 1
    myProj<-paste0("+proj=utm +zone=",myutm," +",hem," +ellps=WGS84")
    myProj
}


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






##' Project altimetry data to the river centerline.
##' access to a basin sqlite file. Basin sqlite files can be acessed from ...  
##' @param datLL a matrix/dataframe with data coordinates (lon,lat)  
##' @param centerlineLL a matrix/dataframe with centerline coordinates (lon,lat)  
##' @importFrom terra vect
##' @importFrom terra project
##' @examples
##' distt<-projectToCenterline(datLL,centerlineLL)
##' out<-getAltReaches(dbfile,mybasin,myreaches)
##' @export

projectToCenterline<-function(datLL,centerLineLL){
   #transform data to UTM coordinates
    myProj<-getProjUTM(centerlineLL[,1], centerlineLL[,2])
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
    return(MyPos)
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
##' @import DBI
##' @import RSQLite
##' @examples
##' dbfile<-"7428.sqlite"
##' mybasin<-'7428'
##' myreaches<-c(74287400183, 74287400163, 74287400151, 74287400173)
##' out<-getAltReaches(dbfile,mybasin,myreaches)
##' @export

getAltReaches<-function(dbfile,mybasin,myreaches){
    mydb <- dbConnect(RSQLite::SQLite(), dbfile)
    out<-lapply(1:length(myreaches),function(i){cat('doing ',i,'\n');getOneReach(mydb,myreaches[i],mybasin)})
    out<-do.call(rbind,out)
    dbDisconnect(mydb)
    out
}










