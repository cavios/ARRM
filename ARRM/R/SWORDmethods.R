# SWORD methods

testCon<-function(tb,ta){
    sqrt((tb$x-ta$x)^2+(tb$y-ta$y)^2)
}

reOrder2<-function(out){
    tt<-split(out,out$reachID)
    tb<-as.numeric(unlist(lapply(tt,nrow)))    
    ta<-rep(1,length(tt))
    for(i in 2:length(tt)){
        nom<-testCon(tt[[i-1]][tb[i-1],],tt[[i]][1,])
        unom<-testCon(tt[[i-1]][tb[i-1],],tt[[i]][tb[i],])
        if(nom > unom) tt[[i]]<-tt[[i]][nrow(tt[[i]]):1,]
    }
    out<-do.call(rbind,tt)
    out
}





reOrder<-function(tt,mytype){
    new<-tt
    if(mytype=='e'){
        new<-tt[nrow(tt):1,]
    }
    new

}



ll2UTM<-function(CLinfo){
    LL<-cbind(CLinfo$x,CLinfo$y)
    myProj<-getProjUTM(LL[,1], LL[,2])
    v <- vect(LL, crs="+proj=longlat")
    u <- terra::project(v, myProj)
    xydat <- crds(u)
    xx<-xydat[,1]
    yy<-xydat[,2]
    cbind(CLinfo,xx,yy)
}


fixCL<-function(CL,lim=300){

    out1<-ll2UTM(CL)
    tt <- split(out1, out1$reachID)
    te <- lapply(tt, getEnd)
    te <-do.call(rbind,te)
    ts <- lapply(tt, getStart)
    ts<-do.call(rbind,ts)
    nr<-length(tt)
    #find end line segment
    firstID<-lapply(1:nrow(te),function(i)findFirst(te[i,],te,ts))
    firstID<-do.call(rbind,firstID)
    idf<-which(firstID$de > lim | firstID$ds > lim)[1]


    # find next line segment
    new<-list()
    myend<-findEnd(te[idf,],te,ts)
    new[[1]]<-tt[[idf]]
    if(myend$dist > lim){
        myend<-findEnd(ts[idf,],te,ts)
        new[[1]]<-reOrder(tt[[idf]],'e')      
    }
    

    nn<-list()
    nn[[1]]<-myend
    for(i in 2:(nr-1)){
        myend<-getNext(myend,te,ts)
        nn[[i]]<-myend
    }
    nn<-do.call(rbind,nn)
    

    
                                        #evt reorder
    id<-nn$id
    for (i in 1:(nr-1)){
        new[[i+1]]<-reOrder(tt[[id[i]]],nn$type[i])      
    }

    
    new<-do.call(rbind,new)
    if(new$wse[1]> new$wse[nrow(new)]) new<-new[nrow(new):1,]
    new
}

findEnd<-function(thiste,te,ts){
    id<-which(te$reachID==thiste$reachID)
    te2<-te[-id,]
    ts2<-ts[-id,]

    ss<-sqrt((thiste$x - ts$x[-id])^2 + (thiste$y - ts$y[-id])^2)
    ee<-sqrt((thiste$x - te$x[-id])^2 + (thiste$y - te$y[-id])^2)
    is<-which.min(ss)
    ie<-which.min(ee)

    if(ss[is] < ee[ie]){
        ids<-which(ts$reachID==ts2[is,1])
        yy<-data.frame(reachid=thiste[1],mymatch=ts2[is,1],id=ids,dist=ss[is],type='s')
    }else{
        ide<-which(te$reachID==te2[ie,1])
        yy<-data.frame(reachid=thiste[1],mymatch=te2[ie,1],id=ide,dist=ee[ie],type='e')
    }
yy

}


findFirst<-function(thiste,te,ts){
    id<-which(te$reachID==thiste$reachID)
    thists<-ts[id,]
    te2<-te[-id,]
    ts2<-ts[-id,]

    ss<-sqrt((thiste$x - ts$x[-id])^2 + (thiste$y - ts$y[-id])^2)
    ee<-sqrt((thiste$x - te$x[-id])^2 + (thiste$y - te$y[-id])^2)
    me<-c(ss,ee)
    ime<-which.min(me)


    ss<-sqrt((thists$x - ts$x[-id])^2 + (thists$y - ts$y[-id])^2)
    ee<-sqrt((thists$x - te$x[-id])^2 + (thists$y - te$y[-id])^2)
    ms<-c(ss,ee)
    ims<-which.min(ms)
    out<-data.frame(reachid=thiste[1],de=me[ime],ds=ms[ims])

    out

}





getNext<-function(myend,te,ts){
    
    if(myend$type=='e'){
        ar1<-ts
    }else{
        ar1<-te
    }
    id<-myend$id
    myend<-findEnd(ar1[id,],te,ts)
    myend
}




getEnd<-function(tt){
    nr<-nrow(tt)
    out<-data.frame(reachID=tt$reachID[nr],x=tt$xx[nr],y=tt$yy[nr])
    out
}

getStart<-function(tt){
    out<-data.frame(reachID=tt$reachID[1],x=tt$xx[1],y=tt$yy[1])
    out
}








getSWORDcode<-function(myarea){
    keySWOT<-data.frame(code=c('1','2','3','4','5','6','7','8','9'),
                        cont=c('af','eu','as','as','oc','sa','na','na','na'))
    idc<-which(keySWOT$code==myarea)
    cont<-keySWOT$cont[idc]
}




##' Extract node information from the SWORD database (Alteneau et al., 2021) for a specific reach
##' Important the function requires the SWORD database! If you do not have this please download from https://www.swordexplorer.com 
##' @param thisreach Reach id from the SWORD database 
##' @param path2sword Path to the location of the SWORD database 
##' @param SWORDv the version of the SWORD database
##' @importFrom hdf5r h5file h5close
##' @examples
##' path2sword<-"/home/karina/SWORD/SWORD_v15_nc/netcdf"
##' SWORDv='15'
##' thisreach<-"21406100031"
##' out<-getReachInfo(thisreachreach,path2sword,SWORDv)
##' @export

getReachInfo<-function(thisreach,path2sword,SWORDv){
    area<-substr(as.character(thisreach),1,1)
    cont<-getSWORDcode(area)     
    file<-h5file(paste0(path2sword,'/',cont,'_sword_v',SWORDv,'.nc'))
    reaches<-file[["nodes/reach_id"]][]
    id<-which(reaches==thisreach)
    out<-NULL
    if(length(id)>0){
        x<-file[["nodes/x"]][id]
        y<-file[["nodes/y"]][id]
        width<-file[["nodes/width"]][id]
        reach_id<-file[["nodes/reach_id"]][id]
        node_id<-file[["nodes/node_id"]][id]
        wse<-file[["nodes/wse"]][id]
        rivername<-file[["nodes/river_name"]][id]
        dist<-file[["nodes/dist_out"]][id]

        out<-data.frame(reachID=reach_id,x=x,y=y,nodeID=node_id,dist=dist,name=rivername,wse=wse)
    }
    h5close(file)
    out
}


getReachInfoDB<-function(thisreach,dbfile,myreach){
    area<-substr(as.character(thisreach),1,1)
    cont<-getSWORDcode(area)     
    dbfile<-paste0(path2DB,'/',cont,'sqlite')
    mydb <- dbConnect(RSQLite::SQLite(), dbfile)
    myq<-paste0("SELECT * FROM nodeinfo"," WHERE reachid =" ,myreach)
    myq<-paste0("SELECT * FROM nodeinfo WHERE reachid = ?", params = list(myreaches))
    out<-dbGetQuery(mydb, "SELECT  * FROM nodeinfo WHERE reachid = ?", params = list(myreaches))

    res <- dbSendQuery(mydb, myq)
    out<-dbFetch(res)
    dbClearResult(res)
    dbDisconnect(mydb)
    out
}


getReachInfo2<-function(dbfile,mybasin,myreaches,addFilter=FALSE,addBuffer=NULL,removeNA=TRUE){
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







##' Extract upstream and downstream reach id from the SWORD database (Alteneau et al., 2021) for a specific reach
##' Important the function requires the SWORD database! If you do not have this please download from https://www.swordexplorer.com 
##' @param thisreach Reach id from the SWORD database 
##' @param path2sword Path to the location of the SWORD database 
##' @param SWORDv the version of the SWORD database
##' @importFrom hdf5r h5file h5close
##' @examples
##' path2sword<-"/home/karina/SWORD/SWORD_v15_nc/netcdf"
##' SWORDv='15'
##' thisreach<-"21406100031"
##' out<-getNeighborReach(thisreachreach,path2sword,SWORDv)
##' @export

getNeighborReach<-function(thisreach,path2sword,SWORDv){
    area<-substr(as.character(thisreach),1,1)
    cont<-getSWORDcode(area)     
    file<-h5file(paste0(path2sword,'/',cont,'_sword_v',SWORDv,'.nc'))
    reaches<-file[["reaches/reach_id"]][]
    id<-which(reaches==thisreach)
    neighborup<-file[["reaches/rch_id_up"]][id,]
    nup<-file[["reaches/n_rch_up"]][id]
    neighbordown<-file[["reaches/rch_id_dn"]][id,]
    out<-c(nup,neighborup,neighbordown)
    h5close(file)
    return(out[out!=0])
}





##' Extract node information from the SWORD database (Alteneau et al., 2021) for multiple reaches and reorders the nodes from downstrean to upstream if these are ordered wrongly
##' Important the function requires the SWORD database! If you do not have this please download from https://www.swordexplorer.com 
##' @param thisreach Reach id from the SWORD database 
##' @param path2sword Path to the location of the SWORD database 
##' @param SWORDv the version of the SWORD database
##' @importFrom hdf5r h5file h5close
##' @examples
##' path2sword<-"/home/karina/SWORD/SWORD_v15_nc/netcdf"
##' SWORDv='15'
##' myreaches<-c(21406100041,21406100031,21406100021)
##' out<-extractCL(myreaches,path2sword,SWORDv)
##' @export

extractCL<-function(myreaches,path2sword,SWORDv){
    o<-order(myreaches)
    out<-lapply(1:length(myreaches),function(i)getReachInfo(myreaches[o[i]],path2sword,SWORDv))
    out<-do.call(rbind,out)
    out<-fixCL(out,lim=300)
    
   out
}


extractCLDB<-function(myreaches,dbfile,lim=300){
    mydb <- dbConnect(RSQLite::SQLite(), dbfile)
    out<-dbGetQuery(mydb, "SELECT  * FROM nodeinfo WHERE reachid = ?", params = list(myreaches))
    dbDisconnect(mydb)
    out<-fixCL(out,lim=lim)
    out
}





