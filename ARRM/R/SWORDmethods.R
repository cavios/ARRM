# SWORD methods

testCon<-function(tb,ta){
    sqrt((tb$x-ta$x)^2+(tb$y-ta$y)^2)
}

reOrder<-function(out){
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

        out<-data.frame(reachID=reach_id,x=x,y=y,nodeID=node_id,dist=dist,name=rivername)
    }
    h5close(file)
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
    o<-order(out$dist)
    out<-out[o,]
    out<-reOrder(out)
    out
}

