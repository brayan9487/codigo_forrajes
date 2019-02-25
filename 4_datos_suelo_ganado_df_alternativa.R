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

# Extraer base de datos de suelo 
soil1 <- raster::stack("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/soil_rasters/TEXMHT_M_sl1_250m.tif")
soil2 <- raster::stack("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/soil_rasters/TEXMHT_M_sl2_250m.tif")
soil3 <- raster::stack("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/soil_rasters/TEXMHT_M_sl3_250m.tif")
soil4 <- raster::stack("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/soil_rasters/TEXMHT_M_sl4_250m.tif")
lista <- list(soil1,soil2,soil3,soil4)

l <- lapply(1:length(lista), function(i){
  cat(paste0("Procesando raster: ",i, "\n" ))
  r  <- lista[[i]] 
  r  <- crop(r,sh)
  r  <- mask(r,sh)
  r  <- raster::resample(r, base)
  df <- rasterToPoints(r)
  df <- as.data.frame(df)
  colnames(df) <- c("lon","lat","Value")
  return(df)
})

l1 <- lapply(2: length(l), function(i){
  x <- l[[i]]
  x<- l[[i]][,-(1:2)]
  return(x)
})
soil  <- do.call(cbind,l1)
soil <- cbind(l[[1]], soil)
rm(soil1,soil2,soil3,soil4)

df <- data.frame(cellID=(cellFromXY(base, xy= soil[,1:2])) , lon= soil$lon ,lat=soil$lat , soil1=soil$Value, soil2= soil$`1`,soil3=soil$`2`,soil4= soil$`3`)

saveRDS(df, "//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/df_soil.rds")

#  Datos dencidad de Ganado, 

cabras   <- raster::stack("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/ganado/cabras/5_Gt_2010_Da.tif")
bufalos  <- raster::stack("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/ganado/bufalos/5_Bf_2010_Da.tif")
ovejas   <- raster::stack("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/ganado/ovejas/5_Sh_2010_Da.tif")
caballos <- raster::stack("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/ganado/caballos/5_Ho_2010_Da.tif")
gana_v   <- raster::stack("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/ganado/ganado_vacuno/5_Ct_2010_Da.tif")
lista <- list(cabras,bufalos,ovejas,caballos,gana_v)
l <- lapply(1:length(lista), function(i){
  r <-lista[[i]]
  r <- crop(r, sh)
  r <- mask(r, sh)
  r <- resample(r, base)
  return(r)
})
r <- raster::stack(l)
df <- as.data.frame(rasterToPoints(r))
df <- data.frame(cellID=(cellFromXY(base,df[,1:2])) , lon=df$x,lat=df$y, cabras=df$X5_Gt_2010_Da, bufalos=df$X5_Bf_2010_Da, ovejas=df$X5_Sh_2010_Da, 
                 caballos=df$X5_Ho_2010_Da, vacuno =df$X5_Ct_2010_Da   )

df <- na.omit(df)

saveRDS(df,"//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/df_ganado.rds")

####
df <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/df_ganado.rds")
test      <- df
test$suma <- NA
test$suma <- test[,4]+test[,5]+test[,6]+test[,7]+test[,8]
test$suma[which(test$suma== 0)] <- NA 
test <- na.omit(test)
df <- test


c <- data.frame(cellID= df$cellID, lon= df$lon,lat=df$lat, cabras=df$cabras)
b <- data.frame(cellID= df$cellID, lon= df$lon,lat=df$lat, bufalos=df$bufalos)
o <- data.frame(cellID= df$cellID, lon= df$lon,lat=df$lat, ovejas=df$ovejas)
ca<- data.frame(cellID= df$cellID, lon= df$lon,lat=df$lat, caballos=df$caballos)
gv<- data.frame(cellID= df$cellID, lon= df$lon,lat=df$lat, gana_v=df$vacuno)


l <- list(c,b,o,ca,gv)
rm(b,c,ca,gv,o,s, samu,test)

samu <- lapply(1:length(l), function(i){
  s <- l[[i]]
  s$condition <- 1
  s$condition[which(s[,4]== 0)]  <- NA
  s <- na.omit(s)
  s$condition <- NULL 
  
  s$condition <- 1
  m <- as.numeric(median(s[,4]))
  s$condition[which(s[,4]< m)] <- NA
  s <- na.omit(s)
  s$condition <- NULL
  colnames(s) <- c("cellID","lon","lat","ganado")
  return(s)
})

bm    <- do.call(rbind,samu)

cellID  <- data.frame(cellID= unique(bm$cellID))
xy <- data.frame(xyFromCell(base, cellID$cellID))

samu <- data.frame(cellID= cellID$cellID, lon= xy$x, lat= xy$y)
saveRDS(samu, "//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/alta_densidad_ganado.rds" )
index <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/vpd_tabla_lista.rds")
rm(cellID,xy)

# que pixeles AÑadir
cellg <- data.frame(cellID= samu$cellID)
celli <- data.frame(cellID= index$cellID)

##
test <- rbind(cellg,celli)
c <- data.frame(unique(test))
nrow(celli)- nrow(cellg)

# resultados 
res <- cellg[which(c(!(cellg$cellID %in% celli$cellID))),]
pix <- data.frame(cellID= res)
xy <- as.data.frame(xyFromCell(base, pix$cellID))

pix <- data.frame(cellID= pix$cellID, lon= xy$x, lat=xy$y)
saveRDS(pix, "//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/agregar_pixeles_alta_densidad_ganado.rds")

### Calcular index para estos pixeles

df   <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/agregar_pixeles_alta_densidad_ganado.rds")
prec <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec.rds")
prec1 <- dplyr::filter(prec, cellID %in% df$cellID)
prec <- prec1

library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)
cores<- detectCores()
cl<- makeCluster(cores-15)
registerDoParallel(cl)

require(parallel)
system.time(indexes_forrages <- foreach(i=1:nrow(prec)) %dopar% {  #nrow(prec_f)
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
  
  
  # 7. P5D: Maximum 5-day running average precipitation (Flooding)
  run_avg <- function(x){
    z <- caTools::runmean(x, k = 5, endrule = 'NA')
    z <- max(z, na.rm = TRUE)
    return(z)
  }
  p5d <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P5D = run_avg(x = Value))
  p5d <- p5d %>% as.data.frame
  names(p5d)[2] <- "Value"; p5d$Variable <- "P5D"
  
  # 8. P_95: 95th percentile of daily precipitation (Erosion risk)
  p_95 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P_95 = quantile(Value, probs = .95, na.rm = TRUE))
  p_95 <- p_95 %>% as.data.frame
  names(p_95)[2] <- "Value"; p_95$Variable <- "P_95"
  
  
  results<- data.frame(cellID= unique(X$cellID), rbind(cdd,totrain,p5d,p_95))
  mylist[[i]] <- results
})

stopCluster(cl)
tabla <- do.call(rbind, indexes_forrages)
saveRDS(tabla,"//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_pixeles_nuevos.rds")

##############################################################################################################

tmin1 <- dplyr::filter(tmin, cellID %in% df$cellID)
tmin <- tmin1

tmax1 <- dplyr::filter(tmax, cellID %in% df$cellID)
tmax <- tmax1

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


######## Index Suplementarios

prec <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/prec_filter_ethiopia.rds")
library(doSNOW)
library(foreach)
library(parallel)
library(doParallel)
cores<- detectCores()
cl<- makeCluster(cores-10)
registerDoParallel(cl)

require(parallel)
system.time(indexes_forrages <- foreach(i=1:nrow(prec)) %dopar% {  #nrow(prec_f)
  require(dplyr)
  mylist <- list()
  time.serie <- prec[i, 1:ncol(prec)]
  
  suppressMessages(library(tidyverse))
  suppressMessages(library(compiler))
  
  X <- time.serie
  X <- X %>% gather(key = Date, value = Value, -(cellID:lat))
  X$Year <- lubridate::year(as.Date(X$Date))
  X$Yday <- lubridate::yday(as.Date(X$Date))
  
  
  # 7. P5D: Maximum 5-day running average precipitation (Flooding)
  run_avg <- function(x){
    z <- caTools::runmean(x, k = 5, endrule = 'NA')
    z <- max(z, na.rm = TRUE)
    return(z)
  }
  p5d <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P5D = run_avg(x = Value))
  p5d <- p5d %>% as.data.frame
  names(p5d)[2] <- "Value"; p5d$Variable <- "P5D"
  
  # 8. P_95: 95th percentile of daily precipitation (Erosion risk)
  p_95 <- X %>% dplyr::group_by(Year) %>% dplyr::arrange(Date) %>% summarise(P_95 = quantile(Value, probs = .95, na.rm = TRUE))
  p_95 <- p_95 %>% as.data.frame
  names(p_95)[2] <- "Value"; p_95$Variable <- "P_95"
  
  
  results<- data.frame(cellID= unique(X$cellID), rbind(p5d,p_95))
  mylist[[i]] <- results
})

stopCluster(cl)
tabla <- do.call(rbind, indexes_forrages)

saveRDS(tabla,"//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_ethiopia_complementarios.rds")

############################################################################################################################################################
colombia <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_colombia_complementarios.rds")
ethiopia <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_ethiopia_complementarios.rds")
kenya    <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_kenya_complementarios.rds")
tanzania <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_tanzania_complementarios.rds")
u_r      <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_uganda_rwanda_complementarios.rds")  

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
saveRDS(tabla2,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_completo_paises_estudio_complementarios.rds")

##### UNion index pasture and crop

df_t <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_completo_paises_estudio.rds")
df_a <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_rain_completo_paises_estudio_complementarios.rds")
index_c <- cbind(df_t,df_a[,-1])
rm(df_t,df_a)

saveRDS(index_c,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cdd_tr_p5_p95_crop_pasture.rds")


##### UNion index ganado
p_a <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_pixeles_nuevos.rds")
unique(p_a$Variable)

vars <- c(as.character(unique(p_a$Variable)))

spread_tables <- lapply(1:length(vars), function(i){
  df <- p_a%>% filter(Variable == vars[i])
  df <- df %>% tidyr::spread(Year, Value) # %>% dplyr::group_by(Variable)
  df$Variable <- NULL
  names(df)[2:ncol(df)] <- paste0(vars[i], "-", names(df)[2:ncol(df)])
  return(df)
})
tabla2 <- Reduce(function(...) merge(..., by = "cellID", all.x= T), spread_tables)


saveRDS(tabla2,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_pixeles_nuevos_completo_paises_estudio.rds")


#########################################################################################################################################################
#### Extraccion de suelo y VPD

library(raster)
c_p <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cdd_tr_p5_p95_crop_pasture.rds")
g   <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/index_pixeles_nuevos_completo_paises_estudio.rds")
index <- rbind(c_p,g)
rm(c_p,g)


soil  <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/df_soil.rds")
xy <- data.frame(xyFromCell(base, index$cellID))

samu <- lapply(1:4, function(i){
  df <- soil 
  df <- df[,-(1:3)]
  df <- df[,i]
  s  <- data.frame(x= soil$lon, y=soil$lat, value=df)
  s  <- rasterFromXYZ(s)
  df1 <- as.data.frame(raster::extract(s,xy))
  return(df1)
})

tabla <- do.call(cbind, samu)
tabla <- cbind(index$cellID, tabla)
colnames(tabla) <- c("cellID","soil1","soil2","soil3","soil4")

index <- cbind(tabla, index[,-1])

saveRDS(index,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_falta_vpd.rds")

###########################################
vpd <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/input_data/5km/vpd_tabla_lista.rds")
xy <- data.frame(xyFromCell(base, index$cellID))

samu <- lapply(1:30, function(i){
  cat(paste0("procesando year",i,"\n"))
  df <- vpd 
  df <- df[,-(1:3)]
  df <- df[,i]
  s  <- data.frame(x= vpd$lon, y=vpd$lat, value=df)
  s  <- rasterFromXYZ(s)
  df1 <- as.data.frame(raster::extract(s,xy))
  return(df1)
})

tabla <- do.call(cbind, samu)
tabla <- cbind(index$cellID, tabla)
colnames(tabla) <- colnames(vpd)

b <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_falta_vpd.rds")
b <- cbind(b, tabla[,-1])
saveRDS(b,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_lista_trabajar.rds")

################################################################################
######################  IDEA DE MEAN, SD, SLOPE  ###############################
################################################################################

##### Definir tabla completa de index para colombia y africa
data_complete  <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_lista_trabajar.rds")  #colombia y este africa

# funcion para calcular pendientes en series de tiempo 
slope <- function(df){
  slope<-numeric(length = nrow(df))
  for(i in 2:(nrow(df))){
    slope[i-1]<-(df[i-1,"val"]-df[i,"val"])/(as.numeric(df[i-1,1]-df[i,1]))
  }
  slope[nrow(df)]<-NA
  df$slope<-slope
  slope<- median(df$slope, na.rm = TRUE) 
  return(slope)
}
# Ejecutamos funcion para reducir tabla completa a mean, sd y slope por index
df <- data_complete
samu <- lapply(1:nrow(df), function(i){
  cat(paste0("Procesando pixel", i, "\n")) 
  year   <- c(1981:2010)
  # suelo 
  cellID <- df[i,1] 
  soil   <- df[i, 2:5] 
  # CDD
  df1         <- t(as.matrix(df[i, 6:35]))
  colnames(df1)  <- c("val")
  rownames(df1)  <- c(1:length(year))
  df1 <- as.data.frame(df1)
  d      <- data.frame( date = year, val = df1$val)
  MCDD      <- mean(d$val)
  SDCDD     <- sd(d$val)
  slopeCDD  <- slope(d)
  #Total Rain 
  df2         <- t(as.matrix(df[i, 36:65]))
  colnames(df2)  <- c("val")
  rownames(df2)  <- c(1:length(year))
  df2 <- as.data.frame(df2)
  d      <- data.frame( date = year, val = df2$val)
  MTR       <- mean(d$val)
  SDTR      <- sd(d$val)
  slopeTR   <- slope(d)
  # P_5D
  df3         <- t(as.matrix(df[i, 66:95]))
  colnames(df3)  <- c("val")
  rownames(df3)  <- c(1:length(year))
  df3 <- as.data.frame(df3)
  d      <- data.frame( date = year, val = df3$val)
  MP_5D       <- mean(d$val)
  SDP_5D      <- sd(d$val)
  slopeP_5D   <- slope(d)
  # P_95D
  df4         <- t(as.matrix(df[i, 96:125]))
  colnames(df4)  <- c("val")
  rownames(df4)  <- c(1:length(year))
  df4 <- as.data.frame(df4)
  d      <- data.frame( date = year, val = df4$val)
  MP_95D       <- mean(d$val)
  SDP_95D      <- sd(d$val)
  slopeP_95D   <- slope(d)
  # VPD 
  df5         <- t(as.matrix(df[i, 126:155]))
  colnames(df5)  <- c("val")
  rownames(df5)  <- c(1:length(year))
  df5 <- as.data.frame(df5)
  d      <- data.frame( date = year, val = df5$val)
  MVPD       <- mean(d$val)
  SDVPD      <- sd(d$val)
  slopeVPD   <- slope(d) 
  return(data.frame(cellID = cellID, soil= soil,
                    mean_CDD= MCDD, sd_CDD= SDCDD, slope_CDD = slopeCDD, 
                    mean_TR= MTR, sd_TR= SDTR, slope_TR= slopeTR,
                    mean_P_5D= MP_5D, sd_P_5D= SDP_5D, slope_P_5D= slopeP_5D,
                    mean_P_95D = MP_95D , sd_P_95D = SDP_95D , slope_P_95D = slopeP_95D ,
                    mean_VPD= MVPD,sd_VPD= SDVPD, slope_VPD= slopeVPD))
})

data_complete <- do.call(rbind,samu)
saveRDS(data_complete, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_completa_mean_sd_slope_lista_trabajar_final.rds")

