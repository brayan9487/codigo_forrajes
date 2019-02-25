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

#Definir Area de estudio
pasture <- raster::stack("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/pasture2000_area.tif")
crop    <- raster::stack("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/cropland2000_area.tif")
pasture <- raster::resample(pasture, base,method= "ngb")
crop    <- raster::resample(crop, base,method= "ngb")
pasture[which(pasture[] <= 0 )] <- NA
pasture[which(pasture[] >   0)] <- 1
crop[which(crop[]<=0)]<- NA
crop[which(crop[]>0)]<-1

crop_area_a <- crop(crop, sh)
crop_area_a <- mask(crop_area_a,sh)

pasture_a <- crop(pasture, sh)
pasture_a <- mask(pasture_a, sh)

crop_area_df <- data.frame(rasterToPoints(crop_area_a))
pasture_df    <- data.frame(rasterToPoints(pasture_a))
rm(crop, pasture, crop_area_a,pasture_a,countries)

crop_area_df$cellID <- cellFromXY(base, xy=data.frame(crop_area_df$x, crop_area_df$y))
pasture_df$cellID <- cellFromXY(base, xy=data.frame(pasture_df$x, pasture_df$y))

names(crop_area_df) <- c("lon","lat","value","cellID")
names(pasture_df) <- c("lon","lat","value","cellID")


cellID1 <- data.frame(crop_area_df$cellID)
names(cellID1) <- "cellID"
cellID2 <- data.frame(pasture_df$cellID)
names(cellID2) <- "cellID"

cellID12 <- rbind (cellID1,cellID2)

area_pasture_crop <- data.frame(cellID= unique(cellID12))
cor<- data.frame(xyFromCell(base, area_pasture_crop$cellID))
rm(cellID1, cellID2 )

area_pasture_crop<- data.frame(cellID= area_pasture_crop$cellID, lon=cor$x, lat = cor$y)
saveRDS(area_pasture_crop, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/crop_pasture.rds")
rm(cor,crop_area_df, pasture_df,cellID12)


#Filtrar bases de datos por crop_area
c_p  <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/crop_pasture.rds") #area de estudio 
prec <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec.rds")
tmax <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmax.rds")
tmin <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmin.rds")

require(dplyr)
prec_f <- dplyr::filter(prec, prec$cellID %in% c_p$cellID)
rm(prec)

tmax_f <- dplyr::filter(tmax, tmax$cellID %in% c_p$cellID)
rm(tmax)

tmin_f <- dplyr::filter(tmin, tmin$cellID %in% c_p$cellID)
rm(tmin)

saveRDS(prec_f, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter.rds")
saveRDS(tmax_f, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmax_filter.rds")
saveRDS(tmin_f, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmin_filter.rds")


# Partir bases por paises
####Correr siempre
countries <- rgdal::readOGR(dsn = paste0(root,"/CWR_pre-breeding/Project_TRUST/Input_data/_shape/_world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
c_a <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/crop_pasture.rds")
prec <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter.rds")

#Colombia
sh_col <- countries[countries@data$COUNTRY== "Colombia",]
b_col <- crop(base, sh_col)
b_col <- mask(b_col, sh_col)
b_col <- data.frame(rasterToPoints(b_col))
b_col <- data.frame(cellID= cellFromXY(base, xy= b_col[,1:2]), lon= b_col$x, lat= b_col$y)
bcol <- dplyr::filter(b_col, b_col$cellID %in% c_a$cellID)
rm(b_col)
prec_col <- dplyr::filter(prec, prec$cellID %in% bcol$cellID)
saveRDS(prec_col, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter_colombia.rds" )

#Ethiopia
sh <- countries[countries@data$COUNTRY== "Ethiopia",]
b_eth <- crop(base, sh)
b_eth <- mask(b_eth,sh)
b_eth <- data.frame(rasterToPoints(b_eth))
b_eth <- data.frame(cellID= cellFromXY(base, xy= b_eth[,1:2]), lon= b_eth$x, lat= b_eth$y)
beth <- dplyr::filter(b_eth, b_eth$cellID %in% c_a$cellID)
rm(b_eth)
prec_eth<- dplyr::filter(prec, prec$cellID %in% beth$cellID)
saveRDS(prec_eth, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter_ethiopia.rds")

#Kenya
sh <- countries[countries@data$COUNTRY== "Kenya",]
b_k <- crop(base, sh)
b_k <- mask(b_k,sh)
b_k <- data.frame(rasterToPoints(b_k))
b_k <- data.frame(cellID= cellFromXY(base, xy= b_k[,1:2]), lon= b_k$x, lat= b_k$y)
bk <- dplyr::filter(b_k, b_k$cellID %in% c_a$cellID)
rm(b_k)
prec_bk<- dplyr::filter(prec, prec$cellID %in% bk$cellID)
saveRDS(prec_bk, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter_kenya.rds")


#Uganda y Ruanda
sh <- countries[countries@data$COUNTRY== "Rwanda"|countries@data$COUNTRY== "Uganda",]
b_ur <- crop(base, sh)
b_ur <- mask(b_ur,sh)
b_ur <- data.frame(rasterToPoints(b_ur))
b_ur <- data.frame(cellID= cellFromXY(base, xy= b_ur[,1:2]), lon= b_ur$x, lat= b_ur$y)
bur <- dplyr::filter(b_ur, b_ur$cellID %in% c_a$cellID)
rm(b_ur)
prec_bur<- dplyr::filter(prec, prec$cellID %in% bur$cellID)
saveRDS(prec_bur, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter_uganda_rwanda.rds")

#Tanzania
sh <- countries[countries@data$COUNTRY== "Tanzania",]
b_t <- crop(base, sh)
b_t <- mask(b_t,sh)

b_t <- data.frame(rasterToPoints(b_t))
b_t <- data.frame(cellID= cellFromXY(base, xy= b_t[,1:2]), lon= b_t$x, lat= b_t$y)

bt <- dplyr::filter(b_t, b_t$cellID %in% c_a$cellID)
rm(b_t)
prec_bt<- dplyr::filter(prec, prec$cellID %in% bt$cellID)
saveRDS(prec_bt, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter_tanzania.rds")



