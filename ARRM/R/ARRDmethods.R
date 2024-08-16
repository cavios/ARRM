

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










