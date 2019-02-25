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

#Calculo de index rain por paises, Windows
prec <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter_tanzania.rds")
library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)
cores<- detectCores()
cl<- makeCluster(cores-31)
registerDoParallel(cl)

require(parallel)
system.time(indexes_forrages <- foreach(i=1:5) %dopar% {  #nrow(prec_f)
  require(dplyr)
  mylist <- list()
  time.serie <- prec[i, 1:ncol(prec)]
  
  suppressMessages(library(tidyverse))
  suppressMessages(library(compiler))
  
  X <- time.serie
  X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
  X$Year <- lubridate::year(as.Date(X$Date))
  X$Yday <- lubridate::yday(as.Date(X$Date))
  
  
  # 1. CDD: Drought spell: Maximum number of consecutive dry days (i.e. with precipitation < 1 mm day-1)
  dr_stress <- function(PREC, p_thresh = 1){
    runs <- rle(PREC < p_thresh)
    cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
    return(cons_days)
  }
  dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
  cdd <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(CDD = dr_stressCMP(Value, p_thresh = 1))
  cdd <- cdd %>% as.data.frame
  names(cdd)[2] <- "Value"; cdd$Variable <- "CDD"
  
  # 2. TOTRAIN: Total precipitation
  X <- time.serie; rm(time.serie)
  X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
  X$Year <- lubridate::year(as.Date(X$Date))
  X$Yday <- lubridate::yday(as.Date(X$Date))
  X$Month <- lubridate::month(as.Date(X$Date))
  
  
  totrain <- X %>% dplyr::group_by(Year)%>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
  totrain <- totrain %>% as.data.frame
  names(totrain)[2] <- "Value"; totrain$Variable <- "TOTRAIN"
  
  
  results<- data.frame(cellID= unique(X$cellID), rbind(cdd,totrain))
  
  mylist[[i]] <- results
})

stopCluster(cl)
tabla <- do.call(rbind, indexes_forrages)
saveRDS(tabla,"//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_uganda_rwanda.rds")


#Calculo de index rain por paises, Linux
prec <- readRDS("/mnt/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter_ethiopia.rds")
require(dplyr)
require(parallel)
system.time(indexes_forrages <- mclapply(1:nrow(prec), function(i){
  cat(paste0("procesando pixel" ,i, "\n"))
  mylist <- list()
  time.serie <- prec[i, 1:ncol(prec)]
  
  suppressMessages(library(tidyverse))
  suppressMessages(library(compiler))
  
  X <- time.serie
  X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
  X$Year <- lubridate::year(as.Date(X$Date))
  X$Yday <- lubridate::yday(as.Date(X$Date))
  
  
  # 1. CDD: Drought spell: Maximum number of consecutive dry days (i.e. with precipitation < 1 mm day-1)
  dr_stress <- function(PREC, p_thresh = 1){
    runs <- rle(PREC < p_thresh)
    cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
    return(cons_days)
  }
  dr_stressCMP <- cmpfun(dr_stress); rm(dr_stress)
  cdd <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(CDD = dr_stressCMP(Value, p_thresh = 1))
  cdd <- cdd %>% as.data.frame
  names(cdd)[2] <- "Value"; cdd$Variable <- "CDD"
  
  # 2. TOTRAIN: Total precipitation
  X <- time.serie; rm(time.serie)
  X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
  X$Year <- lubridate::year(as.Date(X$Date))
  X$Yday <- lubridate::yday(as.Date(X$Date))
  X$Month <- lubridate::month(as.Date(X$Date))
  
  
  totrain <- X %>% dplyr::group_by(Year)%>% dplyr::arrange(Date) %>% summarise(TOTRAIN = sum(Value))
  totrain <- totrain %>% as.data.frame
  names(totrain)[2] <- "Value"; totrain$Variable <- "TOTRAIN"
  results<- data.frame(cellID= unique(X$cellID), rbind(cdd,totrain))
  return(results)
},mc.cores=30, mc.preschedule = F))

saveRDS(tabla,"/mnt/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_ethiopia.rds")


#calculo de VPD 
#Cargar_entradas Linux
tmin <- readRDS("/mnt/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmin.rds")
tmax <- readRDS("/mnt/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/tmax.rds")


require(parallel)
require(dplyr)
system.time(indexes_been <- mclapply(1:nrow(tmax), function(i){   ### nrow(tmax)
  
  cat(paste0("procesando pixel" ,i, "\n"))
  # Parameters
  time.serie  <- tmax[i, 1:ncol(tmax)]
  time.serie1 <- tmin[i, 1:ncol(tmin)]
  
  suppressMessages(library(tidyverse))
  suppressMessages(library(compiler))
  
  X <- time.serie
  X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
  X$Year <- lubridate::year(as.Date(X$Date))
  
  Y <- time.serie1
  Y <- Y %>% gather(key = Date, value = Value, -(cellID:lat))
  Y$Year <- lubridate::year(as.Date(Y$Date))
  
  #Vapour pressure deficit
  calc_vpd <- function(tmin, tmax){
    
    #constants
    albedo <- 0.2
    vpd_cte <- 0.7
    
    #soil heat flux parameters
    a_eslope=611.2
    b_eslope=17.67
    c_eslope=243.5
    
    #input parameters
    tmean <- (tmin+tmax)/2
    
    #soil heat flux
    eslope=a_eslope*b_eslope*c_eslope/(tmean+c_eslope)^2*exp(b_eslope*tmean/(tmean+c_eslope))
    
    #estimate vpd
    esat_min=0.61120*exp((17.67*tmin)/(tmin+243.5))
    esat_max=0.61120*exp((17.67*tmax)/(tmax+243.5))
    vpd=vpd_cte*(esat_max-esat_min) #kPa
    return(vpd)
  }
  
  vpd <- as.data.frame(calc_vpd(tmin = Y$Value, tmax = X$Value ))
  vpd <- data.frame(cellID = X$cellID, lon= X$lon, lat= X$lat, Value= vpd,Year= X$Year) 
  names(vpd)<- c("cellID","lon","lat","Value","Year")
  results <- vpd %>% group_by(Year) %>% summarize(median = median(Value)) 
  results <- data.frame(cellID = unique(vpd$cellID), lon= unique(vpd$lon), lat= unique(vpd$lat),results)
  return(results)
  
}, mc.cores = 40, mc.preschedule = F))

tabla <- do.call(rbind, indexes_been)
saveRDS(taba,"/mnt/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/vpd.rds")

df <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/vpd.rds")
tabla <- df

tabla$Year    <- as.numeric(tabla$Year)
tabla$median   <- as.numeric(tabla$median)
tabla$cellID  <- as.numeric(tabla$cellID)
tabla         <- tabla[complete.cases(tabla),]
tabla         <- unique(tabla)

df <- tabla 
df <- df %>% tidyr::spread(Year, median) 
names(df)[4:ncol(df)] <- paste0("VPD", "-", names(df)[4:ncol(df)])
saveRDS(df,"//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/vpd_tabla_lista.rds")

# Acondicionamientos de tablas de indices 

colombia <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_colombia.rds")
ethiopia <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_ethiopia.rds")
kenya    <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_kenya.rds")
tanzania <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_tanzania.rds")
u_r      <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_uganda_rwanda.rds")  

tabla <- list(colombia,ethiopia,kenya,tanzania,u_r)
rm(colombia,ethiopia,kenya,tanzania,u_r)

tabla2 <- lapply(1:length(tabla), function(j){
  cat(paste0("Processed tabla:", j, "\n"))
  
  tabla[[j]]$Year <- as.numeric(tabla[[j]]$Year)
  tabla[[j]]$Value <- as.numeric(tabla[[j]]$Value)
  tabla[[j]]$cellID <- as.numeric(tabla[[j]]$cellID)
  tabla[[j]] <- tabla[[j]][complete.cases(tabla[[j]]),]
  tabla[[j]] <- unique(tabla[[j]])
  vars <- c(as.character(unique(tabla[[j]]$Variable)))
  
  spread_tables <- lapply(1:length(vars), function(i){
    df <- tabla[[j]] %>% filter(Variable == vars[i])
    df <- df %>% tidyr::spread(Year, Value) # %>% dplyr::group_by(Variable)
    df$Variable <- NULL
    names(df)[2:ncol(df)] <- paste0(vars[i], "-", names(df)[2:ncol(df)])
    return(df)
  })
  tabla2 <- Reduce(function(...) merge(..., by = "cellID", all.x= T), spread_tables)
  return(tabla2)
})
tabla2 <- do.call(rbind, tabla2)
saveRDS(tabla2,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_completo_paises_estudio.rds")