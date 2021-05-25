# carregando pacotes ------------------------------------------------------
library(rgdal);library(LSRS)
library(lattice);library(latticeExtra)
library(dplyr);library(leaflet)
library(gstat);library(sp)
library(leaps);library(Hmisc)
library(raster);library(rgeos)
library(sf);library(units)
library(nngeo);library(mapview)
library(grDevices);
library(RColorBrewer);library(leafsync)
library(reshape);library(corrr)
library(tidyr);library(reshape2)
library(tidyverse);library(caret)
library(glmnet);library(sjPlot)
library(FNN);library(magicfor) 
library(fBasics)

# carregando dados de lab e EM38 e contorno de area 

setwd("D:/OneDrive/Revistas/Soil Systems/salinidade/data/")
crs<-"+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
contorno <- rgdal::readOGR("./contorno_novo.shp")

# plot(contorno)
area_5<-read.asciigrid("./cont_5m.asc",as.image=F,plot.image=T)
area_1<-read.asciigrid("./cont_1m.asc",as.image=F,plot.image=T)
area_10<-read.asciigrid("./boundary_10.asc",as.image=F,plot.image=T)
plot(area_1)
plot(area_5)
plot(area_10)

# Carregando dados de EM38
EM38<-read.table("./teste_em_2.txt", sep=";",header = T,dec = ",")
EM38
# removendo todos os dados negativos de CE.
EM38$CV_1_0m[EM38$CV_1_0m <= 0] <- NA
EM38$CV_0_5m[EM38$CV_0_5m <= 0] <- NA
EM38$IV_1_0m[EM38$IV_1_0m <= 0] <- NA
EM38$IV_0_5m[EM38$IV_0_5m <= 0] <- NA

EM38<-na.omit(EM38)

summary(EM38$CV_1_0m)
summary(EM38$CV_0_5m)
summary(EM38$IV_1_0m)
summary(EM38$IV_0_5m)

# Carregando dados de CE lab
CElab<-read.table("./dados_novos.txt", sep="\t",header = T,dec = ".")
CElab

# Separando dados de CE para validação
CElab_valid <- CElab[c(21:25),]
CElab_train<-dplyr::setdiff(as.data.frame(CElab), CElab_valid)


# estatistica basica 
basicStats(CElab_train[,c(4:6)])
basicStats(CElab_valid[,c(4:6)])
basicStats(CElab[,c(4:6)])
basicStats(EM38[,c(9:12)])

# Histogramas

hist(EM38$CV_1_0m,xlab="aEC 1 m coil",main=NULL)
hist(EM38$CV_0_5m,xlab="aEC 0.5 m coil",main=NULL)
hist(EM38$IV_1_0m,xlab="aMS 1 m coil",main=NULL)
hist(EM38$IV_0_5m,xlab="aMS 0.5 m coil",main=NULL)


hist(log(EM38$CV_1_0m),xlab="log(aEC) 1 m coil",main=NULL)
hist(log(EM38$CV_0_5m),xlab="log(aEC) 0.5 m coil",main=NULL)
hist(log(EM38$IV_1_0m),xlab="log(aMS) 1 m coil",main=NULL)
hist(log(EM38$IV_0_5m),xlab="log(aMS) 0.5 m coil",main=NULL)


# atribuindo coordenadas aos dados de EM38-MK2 e de laboratorio
coordinates(CElab)<-~X+Y
coordinates(CElab_train)<- ~ X+Y
coordinates(CElab_valid)<- ~ X+Y

crs(CElab)<-crs
crs(CElab_train)<-crs
crs(CElab_valid)<-crs
crs(area_1)<-crs
crs(area_5)<-crs
crs(area_10)<-crs

plot(CElab)
plot(CElab_train)
points(CElab_valid,col="red")

# mapa de bolhas 

bub_CE_0_10<-bubble(CElab["CE_0_10"])
bub_CE_10_30<-bubble(CElab["CE_10_30"])
bub_CE_30_50<-bubble(CElab["CE_30_50"])

class(bub_CE_0_10)


bub_CE_lab<-c(bub_CE_0_10,
              bub_CE_10_30,
              bub_CE_30_50,x.same = NA, y.same = NA,
              layout = NULL, merge.legends = FALSE, recursive = FALSE)

bub_CE_lab



# Krigagem dos dados de EM38
coordinates(EM38)<-~X_UTM_24S+Y_UTM_24S
EM38<-EM38[-zerodist(EM38)[,1],]
crs(EM38)<-crs
crs(area_5)<-crs

set.seed(1)
EM_38_size<- floor(0.75 * nrow(EM38)) #25%
set.seed(123)
EM38_retirar <- sample(seq_len(nrow(EM38)), size = EM_38_size)
EM38_2 <- EM38[-EM38_retirar, ] # conjunto de dados - 25%

crs(EM38_2)<-crs

basicStats(data.frame(EM38_2)[,9:12])

par(mfrow=(c(1,2)))
plot(EM38,main="100% de pontos",pch=10,cex=0.2)
plot(EM38_2, main= "25% de pontos",pch=10,cex=0.2)

# CE 1 m
variog.emp.EM_1m<- variogram(log(CV_1_0m) ~ 1, EM38_2, cutoff = 90)
plot(variog.emp.EM_1m)
variog.prelim.EM_1m<- vgm(psill=6200, model="Sph", range=85, nugget=500,boundaries=c(2:3*25, 2:11*50))
variog.aju.EM_1m <- fit.variogram(object=variog.emp.EM_1m, model=variog.prelim.EM_1m, fit.sills=T, fit.ranges=F, fit.method = 7)
plot(variog.emp.EM_1m, variog.aju.EM_1m,main="CEa bobina de 1 m")

# krige.CE_1m<- krigeTg(CV_1_0m~1,locations=EM38_2, newdata=area_10, model=variog.aju.EM_1m,lambda=0)
spplot(krige.CE_1m["var1TG.pred"], main="CEa bobina 1 metro ")
# raster::writeRaster(raster(krige.CE_1m["var1TG.pred"]), filename="D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/CE_1m_TG.tif", overwrite=TRUE,bylayer=T)
CE_1m_TG<-raster("D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/CE_1m_TG.tif")

#CE0.5
variog.emp.EM_0_5m<- variogram(log(CV_0_5m) ~ 1, EM38_2, cutoff = 90)
plot(variog.emp.EM_0_5m)
variog.prelim.EM_0_5m<- vgm(psill=0.35, model="Sph", range=85, nugget=0.07)
variog.aju.EM_0_5m <- fit.variogram(object=variog.emp.EM_0_5m, model=variog.prelim.EM_0_5m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.EM_0_5m, variog.aju.EM_0_5m,main="CEa bobina de 0.5 m")

# krige.CE_0_5m<- krigeTg(CV_0_5m~1,locations=EM38_2, newdata=area_10, model=variog.aju.EM_0_5m,lambda=0)
spplot(krige.CE_0_5m["var1TG.pred"], main="CEa bobina 0.5 metro ")
# raster::writeRaster(raster(krige.CE_0_5m["var1TG.pred"]), filename="D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/CE_0_5m_TG.tif", overwrite=TRUE,bylayer=T)
CE_0_5_m_TG<-raster("D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/CE_0_5m_TG.tif")

#SM1m
variog.emp.SM_1m<- variogram(log(IV_1_0m) ~ 1, EM38_2, cutoff = 180)
plot(variog.emp.SM_1m)
variog.prelim.SM_1m<- vgm(psill=0.15, model="Sph", range=160, nugget=0.005)
variog.aju.SM_1m <- fit.variogram(object=variog.emp.SM_1m, model=variog.prelim.SM_1m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.SM_1m, variog.aju.SM_1m,main="SMa bobina de 1m")

# krige.SM_1m<- krigeTg(IV_1_0m~1,locations=EM38_2, newdata=area_10, model=variog.aju.SM_1m,lambda=0)
spplot(krige.SM_1m["var1TG.pred"], main="SMa bobina 1 metro ")
# raster::writeRaster(raster(krige.SM_1m["var1TG.pred"]), filename="D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/SM_1m_TG.tif", overwrite=TRUE,bylayer=T)
SM_1_m_TG<-raster("D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/SM_1m_TG.tif")

#SM0.5
variog.emp.SM_0_5m<- variogram(log(IV_0_5m) ~ 1, EM38_2, cutoff = 180)
plot(variog.emp.SM_0_5m)
variog.prelim.SM_0_5m<- vgm(psill=0.9, model="Sph", range=150, nugget=0.28)
variog.aju.SM_0_5m <- fit.variogram(object=variog.emp.SM_0_5m, model=variog.prelim.SM_0_5m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.SM_0_5m, variog.aju.SM_0_5m,main="SMa bobina de 0.5 m")

# krige.SM_0_5m<- krigeTg(IV_0_5m~1,locations=EM38_2, newdata=area_10, model=variog.aju.SM_0_5m,lambda=0)
spplot(krige.SM_0_5m["var1TG.pred"], main="SMa bobina 0.5 metro ")
# raster::writeRaster(raster(krige.SM_0_5m["var1TG.pred"]), filename="D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/SM_0_5m_TG.tif", overwrite=TRUE,bylayer=T)
SM_0_5_m_TG<-raster("D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/SM_0_5m_TG.tif")

# plotando semivariogramas juntos 
# tiff(""D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/results/semivariogramas_EM38.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

par(mfrow=c(2,2))

plot(variog.emp.EM_1m$dist, variog.emp.EM_1m$gamma, xlab="Distance (m)",
     ylab="Semivariance",main="aEC (coil 1 m)")

lines(variogramLine(variog.aju.EM_1m, maxdist=100))


plot(variog.emp.EM_0_5m$dist, variog.emp.EM_0_5m$gamma, xlab="Distance (m)",
     ylab="Semivariance",main="aEC (coil 0.5 m)")

lines(variogramLine(variog.aju.EM_0_5m, maxdist=100))

plot(variog.emp.SM_1m$dist, variog.emp.SM_1m$gamma, xlab="Distance (m)",
     ylab="Semivariance",main="aMS (coil 1 m)")

lines(variogramLine(variog.aju.SM_1m, maxdist=200))

plot(variog.emp.SM_0_5m$dist, variog.emp.SM_0_5m$gamma, xlab="Distance (m)",
     ylab="Semivariance",main="aMS (coil 0.5 m)")

lines(variogramLine(variog.aju.SM_0_5m, maxdist=200))

# dev.off()

# Empilhando rasters

sensor<-raster::stack(CE_1m_TG,CE_0_5_m_TG,
                      SM_1_m_TG,SM_0_5_m_TG)

names(sensor)<- c("CE_1_m","CE_0_5_m","SM_1_m","SM_0_5_m")
plot(sensor)

# stackSave(sensor, "D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/stack_sensor_em38.stk")


#Extraindo dados para treinamento
a_sensor_CE_1m = raster::extract(CE_1m_TG, CElab_train, method="simple", sp=T); df_a_sensor_CE_1m = as.data.frame(a_sensor_CE_1m)
a_sensor_CE_05m = raster::extract(CE_0_5_m_TG, CElab_train, method="simple", sp=T); df_a_sensor_CE_05m = as.data.frame(a_sensor_CE_05m)
a_sensor_SM_1m = raster::extract(SM_1_m_TG, CElab_train, method="simple", sp=T); df_a_sensor_SM_1m = as.data.frame(a_sensor_SM_1m)
a_sensor_SM_05m = raster::extract(SM_0_5_m_TG, CElab_train, method="simple", sp=T); df_a_sensor_SM_05m = as.data.frame(a_sensor_SM_05m)

#Extraindo dados para validação
a_sensor_CE_1m_valid = raster::extract(CE_1m_TG, CElab_valid, method="simple", sp=T); df_a_sensor_CE_1m_valid = as.data.frame(a_sensor_CE_1m_valid)
a_sensor_CE_05m_valid = raster::extract(CE_0_5_m_TG, CElab_valid, method="simple", sp=T); df_a_sensor_CE_05m_valid = as.data.frame(a_sensor_CE_05m_valid)
a_sensor_SM_1m_valid = raster::extract(SM_1_m_TG, CElab_valid, method="simple", sp=T); df_a_sensor_SM_1m_valid= as.data.frame(a_sensor_SM_1m_valid)
a_sensor_SM_05m_valid = raster::extract(SM_0_5_m_TG, CElab_valid, method="simple", sp=T); df_a_sensor_SM_05m_valid = as.data.frame(a_sensor_SM_05m_valid)

# Juntando dados de treinamento em tabela 
data_EM_train<-cbind.data.frame(as.data.frame(CElab_train),
                                df_a_sensor_CE_1m[10],df_a_sensor_CE_05m[10],
                                df_a_sensor_SM_1m[10],df_a_sensor_SM_05m[10])
names(data_EM_train)[10:13]<-c("CE_1_m","CE_0_5_m","SM_1_m","SM_0_5_m")

head(data_EM_train)

# write.csv(data_EM_train,"D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/data_ce_lab_pss.csv", row.names = FALSE)

# Juntando dados de validação em tabela 
data_EM_valid<-cbind.data.frame(as.data.frame(CElab_valid),
                                df_a_sensor_CE_1m_valid[10],df_a_sensor_CE_05m_valid[10],
                                df_a_sensor_SM_1m_valid[10],df_a_sensor_SM_05m_valid[10])
names(data_EM_valid)[10:13]<-c("CE_1_m","CE_0_5_m","SM_1_m","SM_0_5_m")

head(data_EM_valid)

# write.csv(data_EM_valid,"D:/OneDrive/EMBRAPA/Pronasolos/Atividades/Publicações/4_diba/diba/data/data_ce_lab_pss_valid.csv", row.names = FALSE)




















































