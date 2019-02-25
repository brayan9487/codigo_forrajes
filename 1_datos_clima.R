#Extraccion de clima 

suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
suppressMessages(if(!require(rgdal)){install.packages('rgdal'); library(rgdal)} else {library(rgdal)})
suppressMessages(if(!require(foreach)){install.packages('foreach'); library(foreach)} else {library(foreach)})
suppressMessages(if(!require(mgcv)){install.packages('mgcv'); library(mgcv)} else {library(mgcv)})
suppressMessages(if(!require(stringr)){install.packages('stringr'); library(stringr)} else {library(stringr)})
suppressMessages(if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)} else {library(tidyverse)})
suppressMessages(if(!require(mapdata)){install.packages('mapdata'); library(mapdata)} else {library(mapdata)})F
suppressMessages(if(!require(xtable)){install.packages('xtable'); library(xtable)} else {library(xtable)})
suppressMessages(if(!require(ggdendro)){install.packages('ggdendro'); library(ggdendro)} else {library(ggdendro)})
suppressMessages(if(!require(compiler)){install.packages('compiler'); library(compiler)} else {library(compiler)})
suppressMessages(if(!require(cluster)){install.packages('cluster'); library(cluster)} else {library(cluster)})
suppressMessages(if(!require(parallelDist)){install.packages('parallelDist'); library(parallelDist)} else {library(parallelDist)})

OSys <- Sys.info(); OSys <- OSys[names(OSys)=="sysname"]
if(OSys == "Linux"){
  root <- "/mnt/workspace_cluster_9"
  base <- raster::stack("/mnt/data_cluster_4/observed/gridded_products/chirps/daily/chirps-v2.0.1981.01.01.tif")
} else {
  if(OSys == "Windows"){
    root <- "//dapadfs/Workspace_cluster_9"
    base <- raster::stack("//dapadfs/data_cluster_4/observed/gridded_products/chirps/daily/chirps-v2.0.1981.01.01.tif")
  }
}; rm(OSys)

countries <- rgdal::readOGR(dsn = paste0(root,"/CWR_pre-breeding/Project_TRUST/Input_data/_shape/_world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sh <- countries[countries@data$COUNTRY== "Colombia"|countries@data$COUNTRY== "Kenya"|countries@data$COUNTRY== "Ethiopia"
                |countries@data$COUNTRY== "Tanzania"|countries@data$COUNTRY== "Uganda"|countries@data$COUNTRY== "Rwanda",]
rm(countries)



# TMIN
library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)

prec <- raster::stack("//dapadfs/data_cluster_5/cropdata/agmerra/daily/nc-files/tmin_daily_ts_agmerra_1980_2010.nc") 
cores<- detectCores()
cl<- makeCluster(cores-7)
registerDoParallel(cl) 

system.time(l <- foreach(i=1:nlayers(prec)) %dopar% {
  mylist <- list()
  
  require(dplyr)
  require(velox)
  require(raster)
  require(raster)
  prec_vx <- velox::velox(prec[[i]])
  pnt <- prec_vx$getCoordinates()
  colnames(pnt) <- c("x", "y")
  pnt <- data.frame(pnt)
  pnt$x[pnt$x > 180] <- pnt$x[pnt$x > 180] - 360
  rtmp <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/base/AgMerra_template.RDS")
  rtmp[which(!is.na(rtmp[]))] <- 1
  pnt_rtmp <- rasterToPoints(rtmp)
  pnt_rtmp <- pnt_rtmp[,1:2]
  pnt_rtmp <- data.frame(pnt_rtmp)
  pnt_rtmp$cellID <- cellFromXY(object = rtmp, xy = pnt_rtmp)
  pnt$cellID <- cellFromXY(object = rtmp, xy = pnt)
  pnt <- pnt[pnt$cellID %in% pnt_rtmp$cellID,]
  rownames(pnt) <- 1:nrow(pnt)
  pnt$x[pnt$x <= 0] <- pnt$x[pnt$x <= 0] + 360
  vls <- prec_vx$extract_points(sp = SpatialPoints(pnt[,1:2]))
  vls <- as_data_frame(vls)
  vls <- cbind(pnt, vls)
  vls$x[vls$x > 180] <- vls$x[vls$x > 180] - 360
  vls1<- data.frame(cellID= vls$cellID, lon= vls$x, lat= vls$y, vls[,-3])
  vls1$x <- NULL 
  vls1$y <- NULL
  r <- rasterFromXYZ(vls1[,-1])
  r <- crop(r, sh) 
  r <- mask(r,sh)
  r <- raster::resample(r, base)
  r  <- rasterToPoints(r)
  r <- data.frame(r)
  r  <- data.frame(cellID= cellFromXY(base, xy= r[,1:2]), lon = r$x, lat= r$y, r$V1)
  mylist[[i]] <- r
} )
stopCluster(cl)

l1 <- lapply(2: length(l), function(i){
  x <- l[[i]]
  x<- l[[i]][,-(1:3)]
  return(x)
})
tmin  <- do.call(cbind,l1)
tmin <- cbind(l[[1]], tmin)
saveRDS(tmin, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmin.rds")

#TMAX
library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)

prec <- raster::stack("//dapadfs/data_cluster_5/cropdata/agmerra/daily/nc-files/tmax_daily_ts_agmerra_1980_2010.nc") 
cores<- detectCores()
cl<- makeCluster(cores-12)
registerDoParallel(cl) 

system.time(l <- foreach(i=1:nlayers(prec)) %dopar% {
  mylist <- list()
  
  require(dplyr)
  require(velox)
  require(raster)
  prec_vx <- velox::velox(prec[[i]])
  pnt <- prec_vx$getCoordinates()
  colnames(pnt) <- c("x", "y")
  pnt <- data.frame(pnt)
  pnt$x[pnt$x > 180] <- pnt$x[pnt$x > 180] - 360
  rtmp <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/base/AgMerra_template.RDS")
  rtmp[which(!is.na(rtmp[]))] <- 1
  pnt_rtmp <- rasterToPoints(rtmp)
  pnt_rtmp <- pnt_rtmp[,1:2]
  pnt_rtmp <- data.frame(pnt_rtmp)
  pnt_rtmp$cellID <- cellFromXY(object = rtmp, xy = pnt_rtmp)
  pnt$cellID <- cellFromXY(object = rtmp, xy = pnt)
  pnt <- pnt[pnt$cellID %in% pnt_rtmp$cellID,]
  rownames(pnt) <- 1:nrow(pnt)
  pnt$x[pnt$x <= 0] <- pnt$x[pnt$x <= 0] + 360
  vls <- prec_vx$extract_points(sp = SpatialPoints(pnt[,1:2]))
  vls <- as_data_frame(vls)
  vls <- cbind(pnt, vls)
  vls$x[vls$x > 180] <- vls$x[vls$x > 180] - 360
  vls1<- data.frame(cellID= vls$cellID, lon= vls$x, lat= vls$y, vls[,-3])
  vls1$x <- NULL 
  vls1$y <- NULL
  r <- rasterFromXYZ(vls1[,-1])
  r <- crop(r, sh) 
  r <- mask(r,sh)
  r <- raster::resample(r, base)
  r  <- rasterToPoints(r)
  r <- data.frame(r)
  r  <- data.frame(cellID= cellFromXY(base, xy= r[,1:2]), lon = r$x, lat= r$y, r$V1)
  mylist[[i]] <- r
} )
stopCluster(cl)
l1 <- lapply(2: length(l), function(i){
  x <- l[[i]]
  x<- l[[i]][,-(1:3)]
  return(x)
})
tmax  <- do.call(cbind,l1)
tmax <- cbind(l[[1]], tmax)
saveRDS(tmax, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmax.rds")

#PREC 
wk_dir <- "//dapadfs/data_cluster_4/observed/gridded_products/chirps/daily"
# Precipitacion
chirps_data <- list.files(wk_dir,pattern='*.tif$', full.names=TRUE) # lista de datos cargados con direccion

library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)

cores<- detectCores()
cl<- makeCluster(cores-5)
registerDoParallel(cl) 

l <- foreach(i=1:length(chirps_data)) %dopar%{
  mylist <- list()
  require(raster)
  r <- raster::stack(chirps_data[[i]])
  r <- crop(r, sh)
  r <- mask(r, sh)
  r <- rasterToPoints(r)
  r <- as.data.frame(r)
  r <- data.frame(cellID= cellFromXY(base, xy= r[,1:2]), lon= r$x, lat= r$y, r[,3])
  mylist[[i]] <- r
}
stopCluster(cl)

l1 <- lapply(2: length(l), function(i){
  x <- l[[i]]
  x<- l[[i]][,-(1:3)]
  return(x)
})
prec  <- do.call(cbind,l1)
prec <- cbind(l[[1]], prec)
saveRDS(prec, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec.rds")

###### Acondicionar tablas 
#prec 
prec_5 <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec.rds")
prec_25 <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/_current_climate/chirps/prec_filtered_oceania.rds")
prec_5 <- prec_5[1:nrow(prec_5), 1:ncol(prec_25)]
colnames(prec_5)<- colnames(prec_25)
prec_5$bean_coordinates <- NULL
saveRDS(prec_5, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec.rds")
rm(prec_5,prec_25)

#Tmax 
tmax_5  <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmax.rds")
tmax_25 <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/_current_climate/agmerra_tmax/tmax_filtered_oceania.rds") 
tmax_5 <- tmax_5[1:nrow(tmax_5), -(4:368)]
colnames(tmax_5) <- colnames(tmax_25)
tmax_5$bean_coordinates <- NULL
saveRDS(tmax_5, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmax.rds")
rm(tmax_25, tmax_5)

#Tmin 
tmin_5  <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmin.rds")
tmin_25 <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/_current_climate/agmerra_tmin/tmin_filtered_oceania.rds")
tmin_5  <- tmin_5[1:nrow(tmin_5), -(4:368)]
colnames(tmin_5) <- colnames(tmin_25)
tmin_5$bean_coordinates <- NULL
saveRDS(tmin_5, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmin.rds")
rm(tmin_25, tmin_5)























