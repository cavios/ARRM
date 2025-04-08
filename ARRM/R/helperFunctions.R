SWORD2RSQLite<-function(cont,path2sword,SWORDv=16){
    file<-h5file(paste0(path2sword,'/',cont,'_sword_v',SWORDv,'.nc'))
    x<-file[["nodes/x"]][]
    y<-file[["nodes/y"]][]
    width<-file[["nodes/width"]][]
    reach_id<-file[["nodes/reach_id"]][]
    node_id<-file[["nodes/node_id"]][]
    wse<-file[["nodes/wse"]][]
    rivername<-file[["nodes/river_name"]][]
    dist<-file[["nodes/dist_out"]][]
    out<-data.frame(reachID=reach_id,nodeID=node_id,x=x,y=y,width=width,wse=wse,dist=dist,rivername=rivername)
    sqlitefile<-paste0(cont,'.sqlite')
    mydb <- dbConnect(RSQLite::SQLite(), sqlitefile)
    tabname<-"nodeinfo"
    dbWriteTable(mydb, tabname, out)
    dbExecute(mydb, 'CREATE INDEX IF NOT EXISTS idx_nodeinfo_reachid ON nodeinfo(reachID)')
    dbDisconnect(mydb)
}


assignPretty <- function(x,n){
  pr <- seq(min(x)-1/365,max(x)+1/365,length=n+1)
  mydiff<-diff(pr)[1]/2
  etatime<-pr[1:n]+mydiff
  br <- pr
  ct <- as.integer(cut(x,br))
  return(list(ct=ct, xx=etatime, br=br))
}


invf <- function(y) -0.5 * log(2/(y + 1) - 1)



createBasinSQliteDB<-function(pathtodat,mybasin){
    cat('reading altimetry data, this may take a while, ..\n')
    dat<-getBasindat(mybasin,pathtodat)
    cat('Creating database, ..\n')
    sqlitefile<-paste0(mybasin,'.sqlite')
    mydb <- dbConnect(RSQLite::SQLite(), sqlitefile)
    tabname<-paste0("basin")
    dbWriteTable(mydb, tabname, dat)
    dbExecute(mydb, 'CREATE INDEX IF NOT EXISTS idx_basin_reachid ON basin(reachID)')
    dbDisconnect(mydb)
}




getAltBasin<-function(myfiles,sat){
    cat('reading data from ',sat,'\n')
    dat<-lapply(1:length(myfiles),function(i)read.table(myfiles[i],header=TRUE))
    dat<-do.call(rbind,dat)
    ids<-which(keySAT$sat==sat)
    dat<-data.frame(timesec=dat$time,time=dat$decyear,lat=dat$lat,lon=dat$lon,height=dat[,keySAT$height[ids]],distnode=dat$distnode, width=dat$width,
                            reachID=dat$reachID, nodeID=dat$nodeID, wse=dat$wse, OCval=dat$OCval,satid=rep(keySAT$satid[ids],nrow(dat)))

    dat
}


getBasindat<-function(mybasin,pathtodat){
    satIC2<-paste0(pathtodat,'/',mybasin,'/IC2')
    if(file.exists(satIC2)){
        myfiles<-dir(paste0(pathtodat,'/',mybasin,'/IC2'),full.names=TRUE)
        ice<-getAltBasin(myfiles,"IC2")
    }else{
        (cat('ICESat-2 data is not available\n'))
    }
    satS3A<-paste0(pathtodat,'/',mybasin,'/S3A')
    
    if(file.exists(satS3A)){
        S3Afiles<-dir(paste0(pathtodat,'/',mybasin,'/S3A'),full.names=TRUE)
        S3A<-getAltBasin(S3Afiles,"S3A")
    }else{
        (cat('S3A data is not available\n'))
    }
    satS3B<-paste0(pathtodat,'/',mybasin,'/S3B')
    if(file.exists(satS3B)){
        S3Bfiles<-dir(paste0(pathtodat,'/',mybasin,'/S3B'),full.names=TRUE)
        S3B<-getAltBasin(S3Bfiles,"S3B")
    }else{
        (cat('S3B data is not available\n'))
    }
    satSAL<-paste0(pathtodat,'/',mybasin,'/SAL')
    
    if(file.exists(satSAL)){
        SALfiles<-dir(paste0(pathtodat,'/',mybasin,'/SAL'),full.names=TRUE)
        SAL<-getAltBasin(SALfiles,"SAL")
    }else{
        (cat('SARAL/AltiKa data is not available\n'))
    }
    satC2<-paste0(pathtodat,'/',mybasin,'/C2')
    if(file.exists(satC2)){
        C2files<-dir(paste0(pathtodat,'/',mybasin,'/C2'),full.names=TRUE)
        C2<-getAltBasin(C2files,"C2")
    }else{
        (cat('Cryosat-2 data is not available\n'))
    }
    
    satS6<-paste0(pathtodat,'/',mybasin,'/S6')
    if(file.exists(satS6)){
        S6files<-dir(paste0(pathtodat,'/',mybasin,'/S6'),full.names=TRUE)
        S6<-getAltBasin(S6files,"S6")
    }else{
        (cat('S6 data is not available\n'))
    }
    datt <- rbind(if(exists("ice")) ice,
                  if(exists("S3A")) S3A,
                  if(exists("S3B")) S3B,
                  if(exists("SAL")) SAL,
                  if(exists("C2")) C2,
                  if(exists("S6")) S6)
    #if(rmNA){
    #    id<-which(!is.na(datt$height))
    #    datt<-datt[id,]
    #}
    datt
}





keySAT<-data.frame(sat=c('C2','SWOT','S3B','S6','SAL','IC2','S3A'),
                   height=c(5,6,6,5,5,5,6),
                   satid=c(1,2,3,4,5,6,7),
                   mycol=c('deeppink1','darkolivegreen3','springgreen2','blue4','deepskyblue1','red3','black')
)




