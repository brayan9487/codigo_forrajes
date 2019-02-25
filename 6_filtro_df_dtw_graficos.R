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

###### AFINAR 
df <- readRDS("//dapadfs/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_lista_trabajar.rds")
df1 <-  df[,36:65]
df <-   cbind(df[,1], df1)
colnames(df) <- c("cellID", colnames(df1))
df1 <- df[,-1]

l <- lapply(1:ncol(df1), function(i){
  cat(paste0("----Procesando columna",i,"\n"))
  d <- df1[,i]
  d[which(d <= 700)] <- 1 
  d <- as.data.frame(d)
  return(d)
})

samu <- do.call(cbind,l)

colnames(samu) <- colnames(df1) 

l <- lapply(1:nrow(samu), function(i){
  cat(paste0("----Procesando pixel",i,"\n"))
  sm <- sum(samu[i,]) 
  return(sm)
})
tr <- as.data.frame(do.call(rbind,l))
samu$condition <- tr
rm(d,df1,tr)
bm <- cbind(df$cellID, samu)
colnames(bm) <- c("cellID", colnames(samu))
bm[which(bm$condition == 30),] <- NA
bm <- na.omit(bm)

cellID_con <- as.data.frame(bm$cellID)
colnames(cellID_con) <- c("cellID")

saveRDS(cellID_con, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cellID_700_mm.rds" )
###########################################################################
africa_cluster   <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cluster_africa_entrenamiento_correccion_4.rds")
colombia_cluster <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cluster_colombia_entrenamiento_correccion_4.rds")
cell_condition   <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cellID_700_mm.rds")


africa_cluster_f    <- dplyr::filter(africa_cluster, cellID %in% cell_condition$cellID)
colombia_cluster_f  <- dplyr::filter(colombia_cluster, cellID %in% cell_condition$cellID)
rm(africa_cluster, colombia_cluster)



#### PLOT Africa
countries <- rgdal::readOGR(dsn = paste0(root,"/CWR_pre-breeding/Project_TRUST/Input_data/_shape/_world_shape"), "all_countries")
proj4string(countries) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

sh_afr <- countries[countries@data$COUNTRY== "Kenya"|countries@data$COUNTRY== "Ethiopia"
                    |countries@data$COUNTRY== "Tanzania"|countries@data$COUNTRY== "Uganda"|countries@data$COUNTRY== "Rwanda",]

cor <- as.data.frame(xyFromCell(base, africa_cluster_f$cellID))
colnames(cor) <- c("lon","lat")
samu <- data.frame(cellID=africa_cluster_f$cellID, lon=cor$lon, lat= cor$lat, cluster= as.factor(as.character(africa_cluster_f$clust)))
samu <- na.omit(samu)

colours <- c("#fdae61", "#d7191c","#2b83ba", "#abdda4")
Y <-ggplot(data = samu, aes(x = lon, y = lat, fill = cluster)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=sh_afr, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.7)+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Environments"))
Y

ggsave(filename="//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cluster_africa_correccion.png",plot=Y, width=10, height=10, units='in')


#### PLOT Colombia
cor <- as.data.frame(xyFromCell(base, colombia_cluster_f$cellID))
colnames(cor) <- c("lon","lat")
samu <- data.frame(cellID=colombia_cluster_f$cellID, lon=cor$lon, lat= cor$lat, cluster= as.factor(as.character(colombia_cluster_f$clust)))
samu <- na.omit(samu)

col <- rgdal::readOGR(dsn ="//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/_shape/col_departamentos_IGAC","Col_dpto_igac_2011_84")
proj4string(col) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


colours <- c("#fdae61", "#d7191c","#2b83ba", "#abdda4")
Y <-ggplot(data = samu, aes(x = lon, y = lat, fill = cluster)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=col, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.7)+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Environments"))
Y

ggsave(filename="//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cluster_colombia_correccion.png",plot=Y, width=10, height=10, units='in')

rm(cell_condition,col,cor,countries,samu,Y)

# pixeles referencia Africa 

a1 <- africa_cluster_f[which(africa_cluster_f$clust== 1),];a1$clust <- NULL ;c1 <- a1$cellID
a2 <- africa_cluster_f[which(africa_cluster_f$clust== 2),];a2$clust <- NULL ;c2 <- a2$cellID
a3 <- africa_cluster_f[which(africa_cluster_f$clust== 3),];a3$clust <- NULL ;c3 <- a3$cellID
a4 <- africa_cluster_f[which(africa_cluster_f$clust== 4),];a4$clust <- NULL ;c4 <- a4$cellID

a1 <- cbind (c1,  a1[,-(1:5)])
a2 <- cbind( c2,  a2[,-(1:5)])
a3 <- cbind( c3,  a3[,-(1:5)])
a4 <- cbind( c4,  a4[,-(1:5)])

r1 <- c(152.83,20.504,0.0000,966.5,153.38, 0.6653,18.000,4.196 ,-0.2045,15.354 ,2.1601 ,-0.01914,1.5358,0.11383,0.010012 )
r2 <- c(71.97 ,18.091, 0.0000 ,1171.9 ,135.00, 1.214 ,19.828 ,3.594, 0.01105 ,18.185 ,2.1388 , 0.06293 ,1.3566  ,0.08602 , 0.004975 )
r3 <- c(165.9 ,32.40,1.000 , 612.7,160.62, -8.813 ,19.26 , 5.436 ,-0.08929,11.293,4.102 , 0.00000 ,1.4817 ,0.08698 ,0.009806  )
r4 <- c(80.47 ,20.866 , 0.0000 ,884.7 ,145.95 , -7.460 ,18.840 , 4.438,-0.04112,14.76,2.5937, 0.0000,1.4238 ,0.08188 , 0.007729  )


euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
dist1 <- lapply(1:nrow(a1) , function(i){
  cat(paste0("procesando pixel ",i,"\n"))
  x <- round(euc.dist(a1[i,-1],r1),3)
  x <- as.data.frame(x)
  names(x) <- c("dist")
  return(x)
})
dist1 <- do.call(rbind, dist1)
head(dist1)


dist2 <- lapply(1:nrow(a2) , function(i){
  cat(paste0("procesando pixel ",i,"\n"))
  x <- round(euc.dist(a2[i,-1],r2),3)
  x <- as.data.frame(x)
  names(x) <- c("dist")
  return(x)
})
dist2 <- do.call(rbind, dist2)
head(dist2)



dist3 <- lapply(1:nrow(a3) , function(i){
  cat(paste0("procesando pixel ",i,"\n"))
  x <- round(euc.dist(a3[i,-1],r3),3)
  x <- as.data.frame(x)
  names(x) <- c("dist")
  return(x)
})
dist3 <- do.call(rbind, dist3)
head(dist3)

dist4 <- lapply(1:nrow(a4) , function(i){
  cat(paste0("procesando pixel ",i,"\n"))
  x <- round(euc.dist(a4[i,-1],r4),3)
  x <- as.data.frame(x)
  names(x) <- c("dist")
  return(x)
})
dist4 <- do.call(rbind, dist4)
head(dist4)
###############################

sort(dist1$dist)
dist1$condition<- NA ;  dist1$condition[which(dist1$dist <= 14.788)] <- 1
dist1<- na.omit(dist1)

sort(dist2$dist)
dist2$condition<- NA ;  dist2$condition[which(dist2$dist <= 15.405)] <- 1
dist2 <- na.omit(dist2)

sort(dist3$dist)
dist3$condition<- NA ;  dist3$condition[which(dist3$dist <= 26.905)] <- 1
dist3 <- na.omit(dist3)

sort(dist4$dist)
dist4$condition<- NA ;  dist4$condition[which(dist4$dist <= 19.601)] <- 1
dist4 <- na.omit(dist4)


r1 <- a1[as.numeric(rownames(dist1)),]
r2 <- a2[as.numeric(rownames(dist2)),]
r3 <- a3[as.numeric(rownames(dist3)),]
r4 <- a4[as.numeric(rownames(dist4)),]


samu <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_lista_trabajar.rds")


ref1 <- dplyr::filter(samu, cellID %in% r1$c1) 
ref2 <- dplyr::filter(samu, cellID %in% r2$c2)
ref3 <- dplyr::filter(samu, cellID %in% r3$c3) 
ref4 <- dplyr::filter(samu, cellID %in% r4$c4) 



d1 = t(apply(ref1[,-(1:5)],2,median,na.rm=TRUE))
p1 <- as.data.frame(d1)
ref1 <- p1


d2 = t(apply(ref2[,-(1:5)],2,median,na.rm=TRUE))
p2 <- as.data.frame(d2)
ref2 <- p2

d3 = t(apply(ref3[,-(1:5)],2,median,na.rm=TRUE))
p3 <- as.data.frame(d3)
ref3 <- p3

d4 = t(apply(ref4[,-(1:5)],2,median,na.rm=TRUE))
p4 <- as.data.frame(d4)
ref4 <- p4

ref1 <- cbind(cellID= 1, ref1)
ref2 <- cbind(cellID= 2, ref2)
ref3 <- cbind(cellID= 3, ref3)
ref4 <- cbind(cellID= 4, ref4)

l <- list(ref1,ref2,ref3,ref4) 
tab_ref <- do.call(rbind,l)

saveRDS(tab_ref, "//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/ref_cellID_clusters_africa.rds")
rm(a1,a2,a3,a4,africa_cluster,cell_condition, cor, d_complet, d1,d2,d3,d4,dist1,dist2,dist3,dist4, l,p1,p2,p3,p4,r1,r2,r3,r4,ref1,ref2,ref3,ref4,samu,x,y)

############# DTW MULTIVARIADA 
samu <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_lista_trabajar.rds")
c1 <- colombia_cluster_f[which(colombia_cluster_f$clust == 1), ]
c2 <- colombia_cluster_f[which(colombia_cluster_f$clust == 2), ]
c3 <- colombia_cluster_f[which(colombia_cluster_f$clust == 3), ]
c4 <- colombia_cluster_f[which(colombia_cluster_f$clust == 4), ]


c1 <- dplyr::filter(samu,cellID %in% c1$cellID)
c2 <- dplyr::filter(samu,cellID %in% c2$cellID)
c3 <- dplyr::filter(samu,cellID %in% c3$cellID)
c4 <- dplyr::filter(samu,cellID %in% c4$cellID)

ref1 <- tab_ref[1,]
ref2 <- tab_ref[2,]
ref3 <- tab_ref[3,]
ref4 <- tab_ref[4,]

##################CLUSTERS Cambiar
cl          <- ref4
table_col   <- c4[,-(2:5)] 

table1 <- rbind(cl,table_col)
#####
tabla2      <- table1
cellID      <- tabla2$cellID

CDD  <- tabla2[,2:31]
CDD  <- cbind(cellID,CDD )

TR  <- tabla2[,32:61]
TR  <- cbind(cellID,TR)

P_5   <- tabla2[,62:91]
P_5   <- cbind(cellID,P_5)

P_95  <- tabla2[,92:121] 
P_95   <- cbind(cellID,P_95)

VPD  <- tabla2[,122:151] 
VPD   <- cbind(cellID,VPD)

require(parallelDist)
dist_mul <- lapply(1:nrow(CDD), function(i){
  cat(paste0("Processed pixel:", i, "\n"))
  
  index_ref_1 <- CDD[1,]
  colnames(index_ref_1) <- as.factor(1:31)
  
  index_ref_2<- TR[1,]
  colnames(index_ref_2) <- as.factor(1:31)
  
  index_ref_3<- P_5[1,]
  colnames(index_ref_3) <- as.factor(1:31)
  
  index_ref_4<- P_95[1,]
  colnames(index_ref_4) <- as.factor(1:31)
  
  index_ref_5<- VPD[1,]
  colnames(index_ref_5) <-as.factor(1:31)
  
  index_ref <- rbind(index_ref_1,index_ref_2,index_ref_3,index_ref_4,index_ref_5)
  
  #### VS 
  index_1 <- CDD[i,]
  colnames(index_1) <- as.factor(1:31)
  
  index_2 <- TR[i,]
  colnames(index_2) <- as.factor(1:31)
  
  index_3 <- P_5[i,]
  colnames(index_3) <-as.factor(1:31)
  
  index_4 <- P_95[i,]
  colnames(index_4) <- as.factor(1:31)
  
  index_5<- VPD[i,]
  colnames(index_5) <- as.factor(1:31)
  
  indexes  <- rbind(index_1,index_2,index_3,index_4,index_5)
  
  df <- rbind(index_ref,indexes)
  
  sample.matrix <-as.matrix(df)
  
  sample.matrix1 <- sample.matrix[1:5,]
  sample.matrix2 <- sample.matrix[6:10,]
  
  sample.matrix1 <- sample.matrix1[,-1]
  sample.matrix2 <- sample.matrix2[,-1]
  sample.List <- list(sample.matrix1, sample.matrix2)
  dist <- parDist(x = sample.List, method = "dtw")
  result <- data.frame(cellID= index_2$`1`, tdist= dist[1] )
  return(result)
})

distance_mult <- do.call(rbind,dist_mul)
###    OJO CAMBIAR 
distance_mult$cluster <- 4
head(distance_mult)
saveRDS(distance_mult,paste0(root,"/CWR_pre-breeding/Project_pastures/Results/distance_dtw_cluster_4.rds"))

#######################################################
colombia_sh <- countries[countries@data$COUNTRY== "Colombia",]

cluster1 <- readRDS(paste0(root,"/CWR_pre-breeding/Project_pastures/Results/distance_dtw_cluster_1.rds"))
cluster2 <- readRDS(paste0(root,"/CWR_pre-breeding/Project_pastures/Results/distance_dtw_cluster_2.rds"))
cluster3 <- readRDS(paste0(root,"/CWR_pre-breeding/Project_pastures/Results/distance_dtw_cluster_3.rds"))
cluster4 <- readRDS(paste0(root,"/CWR_pre-breeding/Project_pastures/Results/distance_dtw_cluster_4.rds"))

tdist <- cluster4
dist2 <-tdist
dist2$cluster <- NULL 
names(dist2) <- c("cellID", "DTWarp")
dist2 <- na.omit(dist2)

dist2$DTWarp_cat <- cut(dist2$DTWarp, quantile(dist2$DTWarp, prob = seq(0,1, by = .05)))
dist2$DTWarp_cat <- as.character(dist2$DTWarp_cat)

library(gtools)
qntls <- as.character(na.omit(unique(dist2$DTWarp_cat)))
qntls <- qntls[gtools::mixedorder(qntls)]
pcnts <- c(paste0(seq(5, 100, 5), "%"))

for(i in 1:length(qntls)){
  dist2$DTWarp_cat[which(dist2$DTWarp_cat == qntls[i])] <- pcnts[i]
}; rm(i)

dist2 <- cbind(dist2, xyFromCell(object = base, cell = dist2$cellID))
dist2$Cat_num <- as.numeric(gsub(pattern = "\\%", replacement = "", dist2$DTWarp_cat))



dist2$DTWarp_cat <- factor(dist2$DTWarp_cat, levels = c("5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"))
dist2$DTWarp_cat <- (as.character(dist2$DTWarp_cat))

##### Puntos de corte de similaridad
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="5%")]  <- "High"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="10%")] <- "High"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="15%")] <- "High"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="20%")] <- "High"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="25%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="30%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="35%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="40%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="45%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="50%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="55%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="60%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="65%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="70%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="75%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="80%")] <- "Medium"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="85%")] <- "Low"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="90%")] <- "Low"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="95%")] <- "Low"
dist2$DTWarp_cat[which(dist2$DTWarp_cat =="100%")] <- "Low"


dist2[which(dist2$DTWarp_cat== "Medium"), ]<- NA
dist2[which(dist2$DTWarp_cat== "Low"), ]<- NA
dist2 <- na.omit(dist2)

saveRDS(dist2, paste0(root,"/CWR_pre-breeding/Project_pastures/Results/5km/high_similarity_cluster4.rds"))


colours <- c("#2ca25f")
col <- rgdal::readOGR(dsn ="//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/_shape/col_departamentos_IGAC","Col_dpto_igac_2011_84")
proj4string(col) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")



cl1  <- readRDS(paste0(root,"/CWR_pre-breeding/Project_pastures/Results/5km/high_similarity_cluster1.rds"))
cl2  <- readRDS(paste0(root,"/CWR_pre-breeding/Project_pastures/Results/5km/high_similarity_cluster2.rds"))
cl3  <- readRDS(paste0(root,"/CWR_pre-breeding/Project_pastures/Results/5km/high_similarity_cluster3.rds"))
cl4  <- readRDS(paste0(root,"/CWR_pre-breeding/Project_pastures/Results/5km/high_similarity_cluster4.rds"))

cl1$con <- "cluster 1" 
cl2$con <- "cluster 2" 
cl3$con <- "cluster 3" 
cl4$con <- "cluster 4" 
l <- list(cl1,cl2 ,cl3 ,cl4 )
df<- do.call(rbind, l)

col <- rgdal::readOGR(dsn ="//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/_shape/col_departamentos_IGAC","Col_dpto_igac_2011_84")
proj4string(col) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

colours <- c("#fdae61", "#d7191c","#2b83ba", "#abdda4")
Y <-ggplot(data = df, aes(x = x, y = y, fill = con)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=col, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.7)+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Environments"))
Y
ggsave(filename="//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/Cluster_colombia_resumen_level_20.png",plot=Y, width=10, height=10, units='in')

