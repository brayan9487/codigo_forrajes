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

d_complet <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/tabla_completa_mean_sd_slope_lista_trabajar_final.rds")

#colombia
sh_col <- countries[countries$COUNTRY== "Colombia",]

r <- crop(base, sh_col); r <- mask(r, sh_col)
r <- as.data.frame(rasterToPoints(r)); r <-r[,-3]
colnames(r) <-c("lon","lat")
r<- data.frame(cellID= (cellFromXY(base, r)), lon= r$lon, lat=r$lat)
colombia <- dplyr::filter(d_complet, cellID %in% r$cellID);rm(r,sh_col)

#extraer africa
sh_afr <- countries[countries@data$COUNTRY== "Kenya"|countries@data$COUNTRY== "Ethiopia"
                    |countries@data$COUNTRY== "Tanzania"|countries@data$COUNTRY== "Uganda"|countries@data$COUNTRY== "Rwanda",]
r <- crop(base, sh_afr); r <- mask(r, sh_afr)
r <- as.data.frame(rasterToPoints(r)); r <-r[,-3]
colnames(r) <-c("lon","lat")
r<- data.frame(cellID= (cellFromXY(base, r)), lon= r$lon, lat=r$lat)
africa <- dplyr::filter(d_complet, cellID %in% r$cellID);rm(r)
africa <- na.omit(africa)


#muestras
df  <- africa[,-1]
df1 <- round(df[,1:4],digits=0)
df <-  cbind(df1, df[,5:ncol(df)]);rm(df1)

df1 <-as.data.frame(as.factor(df[,1]))
df2 <-as.data.frame(as.factor(df[,2]))
df3 <-as.data.frame(as.factor(df[,3]))
df4 <-as.data.frame(as.factor(df[,4]))

d <- cbind(df1,df2,df3,df4, df[,5:ncol(df)])
colnames(d) <- c("soil1","soil2","soil3","soil4", colnames(df[,5:ncol(df)]))
str(d)
rm(df1,df2,df3,df4,df)

d <- cbind(africa$cellID, d)
colnames(d) <- c("cellID", colnames(d[,-1]))

set.seed(123)
s_africa <-  d[sample(1:nrow(d),5000,replace = FALSE),]
s_cell <- s_africa$cellID
import <- s_africa

## Matriz de distancias que necesito 
system.time( daisy.mat <- as.matrix(daisy(s_africa[,-1], metric="gower", type=list(ordratio=c(1,2,3,4)))))
a <- as.dist(daisy.mat)
clust_hc <- fastcluster::hclust(a, method = "ward.D"); rm(daisy.mat)
hcd <- as.dendrogram(clust_hc)
plot(cut(hcd, h = 10)$upper, main="Upper tree of cut at h=75")
# abline(h=150,col="red")
memb <- cutree(clust_hc, k= 4)


d_clust <- data.frame( s_africa , clust= factor(memb) )
d_clust$clust <- as.factor(d_clust$clust)
rm(d,daisy.mat, df,d)


africa_entrenamiento <- d_clust 

library(randomForest)
system.time(Modelo<-randomForest(clust ~ ., 
                                 data=africa_entrenamiento[,-1], # datos para entrenar 
                                 ntree=1000,           # cantidad de arboles   
                                 mtry=6,             # cantidad de variables
                                 replace=T)  )        # muestras con reemplazo

s_cell <- data.frame(cellID= s_cell)
africa_prueba <- africa[!(africa$cellID %in% s_cell$cellID), ]
cell_prueba <- africa_prueba$cellID
df  <- africa_prueba[,-1]
df1 <- round(df[,1:4],digits=0)
df <-  cbind(df1, df[,5:ncol(df)]);rm(df1)

df1 <-as.data.frame(as.factor(df[,1]))
df2 <-as.data.frame(as.factor(df[,2]))
df3 <-as.data.frame(as.factor(df[,3]))
df4 <-as.data.frame(as.factor(df[,4]))

d <- cbind(df1,df2,df3,df4, df[,5:ncol(df)])
colnames(d) <- c("soil1","soil2","soil3","soil4", colnames(df[,5:ncol(df)]))
str(d)
rm(df1,df2,df3,df4,df)
africa_prueba <- d 
str(africa_prueba)
# levels(africa_prueba$soil3) <- c(1:10)

Prediccion <- predict(Modelo , africa_prueba)
table(Prediccion)

africa_prueba <- cbind(cell_prueba, africa_prueba, Prediccion)

colnames(africa_prueba) <- colnames(africa_entrenamiento)
africa_cluster <- rbind(africa_entrenamiento,africa_prueba)
saveRDS(africa_cluster,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cluster_africa_entrenamiento_correccion_4.rds")

#
cell_colombia <- colombia$cellID
df  <- colombia[,-1]
df1 <- round(df[,1:4],digits=0)
df <-  cbind(df1, df[,5:ncol(df)]);rm(df1)

df1 <-as.data.frame(as.factor(df[,1]))
df2 <-as.data.frame(as.factor(df[,2]))
df3 <-as.data.frame(as.factor(df[,3]))
df4 <-as.data.frame(as.factor(df[,4]))

d <- cbind(df1,df2,df3,df4, df[,5:ncol(df)])
colnames(d) <- c("soil1","soil2","soil3","soil4", colnames(df[,5:ncol(df)]))
str(d)
rm(df1,df2,df3,df4,df)
levels(d$soil1) <- c(1:10)
levels(d$soil2) <- c(1:10)
levels(d$soil3) <- c(1:10)

Prediccion <- predict (Modelo , d)
table(Prediccion)

colombia_cluster <- cbind(cell_colombia, d, Prediccion)
colnames(colombia_cluster)  <- colnames(africa_cluster)
saveRDS(colombia_cluster,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cluster_colombia_entrenamiento_correccion_4.rds")
# 
rm(africa, africa_entrenamiento, africa_prueba, colombia, d,d_clust, d_complet, hcd, s_africa, s_cell
)


#### PLOT Africa
sh_afr <- countries[countries@data$COUNTRY== "Kenya"|countries@data$COUNTRY== "Ethiopia"
                    |countries@data$COUNTRY== "Tanzania"|countries@data$COUNTRY== "Uganda"|countries@data$COUNTRY== "Rwanda",]


africa_cluster <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cluster_africa_entrenamiento_correccion_4.rds")
cor <- as.data.frame(xyFromCell(base, africa_cluster$cellID))
colnames(cor) <- c("lon","lat")
samu <- data.frame(cellID=africa_cluster$cellID, lon=cor$lon, lat= cor$lat, cluster= as.factor(as.character(africa_cluster$clust)))
samu <- na.omit(samu)

colours <- c("#d7191c", "#fdae61","#abdda4","#2b83ba")
Y <-ggplot(data = samu, aes(x = lon, y = lat, fill = cluster)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=sh_afr, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.1)+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Environments"))
Y

#### PLOT Colombia

colombia_cluster <- readRDS("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Project_pastures/results/5km/cluster_colombia_entrenamiento_correccion_4.rds")

cor <- as.data.frame(xyFromCell(base, colombia_cluster$cellID))
colnames(cor) <- c("lon","lat")
samu <- data.frame(cellID=colombia_cluster$cellID, lon=cor$lon, lat= cor$lat, cluster= as.factor(as.character(colombia_cluster$clust)))
samu <- na.omit(samu)

col <- rgdal::readOGR(dsn ="//dapadfs.cgiarad.org/workspace_cluster_9/CWR_pre-breeding/Project_TRUST/Input_data/_shape/col_departamentos_IGAC","Col_dpto_igac_2011_84")
proj4string(col) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


colours <- c("#d7191c", "#fdae61","#abdda4","#2b83ba")
Y <-ggplot(data = samu, aes(x = lon, y = lat, fill = cluster)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual (values=colours,na.value = "gray") + 
  geom_polygon(data=col, aes(x=long, y=lat, group=group),fill=NA,color="black", size=0.1)+
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill = guide_legend(title = "Environments"))
Y

rm(cor,col,samu) 
