# carregando pacotes ------------------------------------------------------
library(rgdal)
library(LSRS)
library(lattice)
library(latticeExtra)
library(dplyr)
library(leaflet)
library(gstat)
library(sp)
library(leaps)
library(Hmisc)
library(raster)
library(rgeos)
library(sf)
library(units)
library(nngeo)
library(mapview)
library(grDevices)
library(RColorBrewer)
library(leafsync)
library(reshape)
library(corrr)
library(tidyr)
library(reshape2)
library(tidyverse)
library(caret)
library(glmnet)
library(sjPlot)
library(FNN)
library(magicfor) 
library(fBasics)


# carregando dados de lab e EM38 e contorno de area ------------------------------------------------------

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

CElab<-read.table("./dados_novos.txt", sep="\t",header = T,dec = ".")
CElab

CElab_valid <- CElab[c(21:25),]
CElab_train<-dplyr::setdiff(as.data.frame(CElab), CElab_valid)

# CElab[, c(2:9)] <- sapply(CElab[, c(2:9)], as.numeric)
# CElab_train[, c(2:7)] <- sapply(CElab_train[, c(2:7)], as.numeric)
# CElab_valid[, c(2:7)] <- sapply(CElab_valid[, c(2:7)], as.numeric)

# estatistica basica ------------------------------------------------------
basicStats(CElab_train[,c(4:6)])
basicStats(CElab_valid[,c(4:6)])
basicStats(CElab[,c(4:6)])


basicStats(EM38[,c(9:12)])



# Histogramas -------------------------------------------------------------

hist(EM38$CV_1_0m,xlab="aEC 1 m coil",main=NULL)
hist(EM38$CV_0_5m,xlab="aEC 0.5 m coil",main=NULL)
hist(EM38$IV_1_0m,xlab="aMS 1 m coil",main=NULL)
hist(EM38$IV_0_5m,xlab="aMS 0.5 m coil",main=NULL)


hist(log(EM38$CV_1_0m),xlab="log(aEC) 1 m coil",main=NULL)
hist(log(EM38$CV_0_5m),xlab="log(aEC) 0.5 m coil",main=NULL)
hist(log(EM38$IV_1_0m),xlab="log(aMS) 1 m coil",main=NULL)
hist(log(EM38$IV_0_5m),xlab="log(aMS) 0.5 m coil",main=NULL)


# atribuindo coordenadas aos dados de EM38-MK2 e de laboratorio -----------
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

# mapa de bolhas ----------------------------------------------------------

bub_CE_0_10<-bubble(CElab["CE_0_10"])
bub_CE_10_30<-bubble(CElab["CE_10_30"])
bub_CE_30_50<-bubble(CElab["CE_30_50"])

class(bub_CE_0_10)


bub_CE_lab<-c(bub_CE_0_10,
              bub_CE_10_30,
              bub_CE_30_50,x.same = NA, y.same = NA,
              layout = NULL, merge.legends = FALSE, recursive = FALSE)

bub_CE_lab



# Krig EM38 ---------------------------------------------------------------
coordinates(EM38)<-~X_UTM_24S+Y_UTM_24S
EM38<-EM38[-zerodist(EM38)[,1],]
crs(EM38)<-crs
crs(area_5)<-crs

set.seed(1)
EM_38_size<- floor(0.95 * nrow(EM38)) #25%
set.seed(123)
EM38_retirar <- sample(seq_len(nrow(EM38)), size = EM_38_size)
EM38_2 <- EM38[-EM38_retirar, ] # conjunto de dados - 25%

crs(EM38_2)<-crs

basicStats(data.frame(EM38_2)[,9:12])

par(mfrow=(c(1,2)))
plot(EM38,main="100% ou 20669 pontos",pch=10,cex=0.2)
plot(EM38_2, main= "5% ou 1034 pontos",pch=10,cex=0.2)

# CE 1 m
variog.emp.EM_1m<- variogram(log(CV_1_0m) ~ 1, EM38_2, cutoff = 90)
plot(variog.emp.EM_1m)
variog.prelim.EM_1m<- vgm(psill=6200, model="Sph", range=85, nugget=500,boundaries=c(2:3*25, 2:11*50))
variog.aju.EM_1m <- fit.variogram(object=variog.emp.EM_1m, model=variog.prelim.EM_1m, fit.sills=T, fit.ranges=F, fit.method = 7)
plot(variog.emp.EM_1m, variog.aju.EM_1m,main="CEa bobina de 1 m")

# krige.CE_1m<- krigeTg(CV_1_0m~1,locations=EM38_2, newdata=area_5, model=variog.aju.EM_1m,lambda=0)
# spplot(krige.CE_1m["var1TG.pred"], main="CEa bobina 1 metro ")

# raster::writeRaster(raster(krige.CE_1m["var1TG.pred"]), filename="D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_1m_TG.tif", overwrite=TRUE,bylayer=T)
CE_1m_TG<-raster("D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_1m_TG.tif")


#CE0.5
variog.emp.EM_0_5m<- variogram(log(CV_0_5m) ~ 1, EM38_2, cutoff = 90)
plot(variog.emp.EM_0_5m)
variog.prelim.EM_0_5m<- vgm(psill=0.42, model="Sph", range=85, nugget=0.05)
variog.aju.EM_0_5m <- fit.variogram(object=variog.emp.EM_0_5m, model=variog.prelim.EM_0_5m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.EM_0_5m, variog.aju.EM_0_5m,main="CEa bobina de 0.5 m")

# krige.CE_0_5m<- krigeTg(CV_0_5m~1,locations=EM38_2, newdata=area_5, model=variog.aju.EM_0_5m,lambda=0)
# spplot(krige.CE_0_5m["var1TG.pred"], main="CEa bobina 0.5 metro ")

# raster::writeRaster(raster(krige.CE_0_5m["var1TG.pred"]), filename="D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_0_5m_TG.tif", overwrite=TRUE,bylayer=T)
CE_0_5_m_TG<-raster("D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_0_5m_TG.tif")


#SM1m
variog.emp.SM_1m<- variogram(log(IV_1_0m) ~ 1, EM38_2, cutoff = 180)
plot(variog.emp.SM_1m)
variog.prelim.SM_1m<- vgm(psill=0.15, model="Sph", range=160, nugget=0.005)
variog.aju.SM_1m <- fit.variogram(object=variog.emp.SM_1m, model=variog.prelim.SM_1m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.SM_1m, variog.aju.SM_1m,main="SMa bobina de 1m")

# krige.SM_1m<- krigeTg(IV_1_0m~1,locations=EM38_2, newdata=area_5, model=variog.aju.SM_1m,lambda=0)
# spplot(krige.SM_1m["var1TG.pred"], main="SMa bobina 1 metro ")

# raster::writeRaster(raster(krige.SM_1m["var1TG.pred"]), filename="D:/OneDrive/Revistas/Soil Systems/salinidade/data/SM_1m_TG.tif", overwrite=TRUE,bylayer=T)
SM_1_m_TG<-raster("D:/OneDrive/Revistas/Soil Systems/salinidade/data/SM_1m_TG.tif")


#SM0.5
variog.emp.SM_0_5m<- variogram(log(IV_0_5m) ~ 1, EM38_2, cutoff = 180)
plot(variog.emp.SM_0_5m)
variog.prelim.SM_0_5m<- vgm(psill=0.9, model="Sph", range=150, nugget=0.28)
variog.aju.SM_0_5m <- fit.variogram(object=variog.emp.SM_0_5m, model=variog.prelim.SM_0_5m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.SM_0_5m, variog.aju.SM_0_5m,main="SMa bobina de 0.5 m")

# krige.SM_0_5m<- krigeTg(IV_0_5m~1,locations=EM38_2, newdata=area_5, model=variog.aju.SM_0_5m,lambda=0)
# spplot(krige.SM_0_5m["var1TG.pred"], main="SMa bobina 0.5 metro ")

# raster::writeRaster(raster(krige.SM_0_5m["var1TG.pred"]), filename="D:/OneDrive/Revistas/Soil Systems/salinidade/data/SM_0_5m_TG.tif", overwrite=TRUE,bylayer=T)
SM_0_5_m_TG<-raster("D:/OneDrive/Revistas/Soil Systems/salinidade/data/SM_0_5m_TG.tif")

# plotando semivariogramas juntos -----------------------------------------
tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/semivariogramas_EM38.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

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

dev.off()


# Empilhando rasters ------------------------------------------------------

sensor<-raster::stack(CE_1m_TG,CE_0_5_m_TG,
                      SM_1_m_TG,SM_0_5_m_TG)

names(sensor)<- c("CE_1_m","CE_0_5_m","SM_1_m","SM_0_5_m")


tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/CE_1m.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

plot(sensor[[1]],col=hcl.colors(20, palette = "Earth",rev=T))

dev.off()

tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/CE_05m.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

plot(sensor[[2]],col=hcl.colors(20, palette = "Earth",rev=T))

dev.off()

tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/SM_1m.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

plot(sensor[[3]],col=hcl.colors(20, palette = "Earth",rev=T))

dev.off()

tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/SM_05m.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

plot(sensor[[4]],col=hcl.colors(20, palette = "Earth",rev=T))

dev.off()

tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/EM38_juntos.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

plot(sensor,col=hcl.colors(20, palette = "Earth",rev=T),main=c("aEC (1 m coil)","aEC (0.5 m coil)","aMS (1 m coil)","aMS (0.5 m coil)"))

dev.off()

points <- SpatialPoints(coords = CElab_train, proj4string = CRS("+proj=longlat +datum=WGS84"))
points_valid <- SpatialPoints(coords = CElab_valid, proj4string = CRS("+proj=longlat +datum=WGS84"))

coordinates(CElab)<-~X+Y

points <- SpatialPoints(coords = CElab, proj4string = CRS("+proj=longlat +datum=WGS84"))

fun <- function() {
  plot(points, add = TRUE, col = "red", pch = 3)
  plot(points_valid, add = TRUE, col = "blue", pch = 3)
}

fun <- function() {
  plot(points, add = TRUE, col = "red", pch = 3)
}


tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/EM38_juntos_pontos_CE_lab.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

plot(sensor,addfun=fun,col=hcl.colors(20, palette = "Earth",rev=T),main=c("aEC (1 m coil)","aEC (0.5 m coil)","aMS (1 m coil)","aMS (0.5 m coil)"))

dev.off()

sensor_grid <- as(sensor, "SpatialPixelsDataFrame")
sensor_grid <- na.omit(sensor_grid)

# extrair dados para coordenadas e empilhando ------------------------------------------
#treinamento

a_sensor_CE_1m = raster::extract(CE_1m_TG, CElab, method="simple", sp=T); df_a_sensor_CE_1m = as.data.frame(a_sensor_CE_1m)
a_sensor_CE_05m = raster::extract(CE_0_5_m_TG, CElab, method="simple", sp=T); df_a_sensor_CE_05m = as.data.frame(a_sensor_CE_05m)
a_sensor_SM_1m = raster::extract(SM_1_m_TG, CElab, method="simple", sp=T); df_a_sensor_SM_1m = as.data.frame(a_sensor_SM_1m)
a_sensor_SM_05m = raster::extract(SM_0_5_m_TG, CElab, method="simple", sp=T); df_a_sensor_SM_05m = as.data.frame(a_sensor_SM_05m)


data_EM_train<-cbind.data.frame(as.data.frame(CElab),
                                df_a_sensor_CE_1m[10],df_a_sensor_CE_05m[10],
                                df_a_sensor_SM_1m[10],df_a_sensor_SM_05m[10])
names(data_EM_train)[10:13]<-c("CE_1_m","CE_0_5_m","SM_1_m","SM_0_5_m")

head(data_EM_train)



# correlacao --------------------------------------------------------------
cor_EM_train<-cor(data_EM_train[,c(4:6,10:13)],use="pairwise.complete.obs")

cor_EM_train

# data_EM_train[, c(2:13)] <- sapply(data_EM_train[, c(2:13)], as.numeric)
# data_EM_valid[, c(2:13)] <- sapply(data_EM_valid[, c(2:13)], as.numeric)

# Regressoes --------------------------------------------------------------
#CE 0-10
EC_0_10<-regsubsets(CE_0_10 ~ .,data=data_EM_train[,c(4,10:13)],nbest = 50, nvmax = 50,method = "exhaustive",really.big=T)

EC_0_10_reg<- summary(EC_0_10)

coef(EC_0_10, id=which(EC_0_10_reg$bic == min(EC_0_10_reg$bic)))
coef(EC_0_10, id=which(EC_0_10_reg$adjr2 == max(EC_0_10_reg$adjr2)))

EC_0_10_lm_r2<-lm(CE_0_10 ~ CE_0_5_m,data=data_EM_train)
summary(EC_0_10_lm_r2)

#CE 10-30
EC_10_30<-regsubsets(CE_10_30 ~ .,data=data_EM_train[,c(5,10:13)],nbest = 50, nvmax = 50,method = "exhaustive",really.big=T)

EC_10_30_reg<- summary(EC_10_30)

coef(EC_10_30, id=which(EC_10_30_reg$bic == min(EC_10_30_reg$bic)))
coef(EC_10_30, id=which(EC_10_30_reg$adjr2 == max(EC_10_30_reg$adjr2)))

EC_10_30_lm_r2<-lm(CE_10_30 ~ CE_0_5_m+SM_1_m+SM_0_5_m ,data=data_EM_train)
summary(EC_10_30_lm_r2)

#CE 10-30
EC_30_50<-regsubsets(CE_30_50 ~ .,data=data_EM_train[,c(6,10:13)],nbest = 50, nvmax = 50,method = "exhaustive",really.big=T)
EC_30_50_reg<- summary(EC_30_50)

coef(EC_30_50, id=which(EC_30_50_reg$bic == min(EC_30_50_reg$bic)))
coef(EC_30_50, id=which(EC_30_50_reg$adjr2 == max(EC_30_50_reg$adjr2)))

EC_30_50_lm_r2<-lm(CE_30_50 ~  CE_1_m+CE_0_5_m+SM_1_m,data=data_EM_train)
summary(EC_30_50_lm_r2)

# plotando modelos juntos -------------------------------------------------

tab_model(
  EC_0_10_lm_r2,EC_10_30_lm_r2,EC_30_50_lm_r2,
  dv.labels = c("CE (pred 0-10 cm) ", "CE (pred 10-30 cm)","CE (pred 30-50 cm)"),
  string.pred = "Coeffcient",
  string.ci = "Conf. Int (95%)",
  string.p = "P-Value"
)


# predict -----------------------------------------------------------------
CE_0_10_pred<-raster::predict(sensor,EC_0_10_lm_r2)
crs(CE_0_10_pred)<-"+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
names(CE_0_10_pred)<-"CE_0_10_reg"
CE_0_10_pred[CE_0_10_pred$CE_0_10_reg <= 0]<-0
spplot(CE_0_10_pred) 

CE_10_30_pred<-predict(sensor,EC_10_30_lm_r2)
crs(CE_10_30_pred)<-"+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
names(CE_10_30_pred)<-"CE_10_30_reg"
CE_10_30_pred[CE_10_30_pred$CE_10_30_reg <= 0]<-0
spplot(CE_10_30_pred)

CE_30_50_pred<-predict(sensor,EC_30_50_lm_r2)
crs(CE_30_50_pred)<-"+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
names(CE_30_50_pred)<-"CE_30_50_reg"
CE_30_50_pred[CE_30_50_pred$CE_30_50_reg <= 0]<-0
spplot(CE_30_50_pred)

CE_pred<-addLayer(CE_0_10_pred,
                  CE_10_30_pred,
                  CE_30_50_pred)
tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/EC_pred_abordagem1.tiff", height = 10, width = 18, units="cm",compression = "lzw", res = 300)

spplot(CE_pred,col.regions=hcl.colors(50, palette = "RdYlBu",rev=T),
       names.attr=c('(pred) EC lab 0-10 cm','(pred) EC lab10-30 cm','(pred) EC lab 30-50 cm'),
       at= seq(0, 45, 1))

dev.off()

plot(sensor)



# preparacao para validacao -----------------------------------------------

a_CE_0_10_cm = raster::extract(CE_0_10_pred, CElab, method="simple", sp=T); df_CE_0_10_cm = as.data.frame(a_CE_0_10_cm)
a_CE_10_30_cm = raster::extract(CE_10_30_pred, CElab, method="simple", sp=T); df_CE_10_30_cm = as.data.frame(a_CE_10_30_cm)
a_CE_30_50_cm = raster::extract(CE_30_50_pred, CElab, method="simple", sp=T); df_CE_30_50_cm = as.data.frame(a_CE_30_50_cm)

# validacao ---------------------------------------------------------------

res_CE_0_10<-mean(df_CE_0_10_cm$CE_0_10 - df_CE_0_10_cm$CE_0_10_reg)
res_CE_0_10

res_CE_10_30<-mean(df_CE_10_30_cm$CE_10_30 - df_CE_10_30_cm$CE_10_30_reg)
res_CE_10_30

res_CE_30_50<-mean(df_CE_30_50_cm$CE_30_50 - df_CE_30_50_cm$CE_30_50_reg)
res_CE_30_50

rmse_CE_0_10 <- sqrt(sum((df_CE_0_10_cm$CE_0_10_reg - df_CE_0_10_cm$CE_0_10)^2) / nrow(df_CE_0_10_cm))
rmse_CE_0_10 #quanto menor, melhor

rmse_CE_10_30 <- sqrt(sum((df_CE_10_30_cm$CE_10_30_reg - df_CE_10_30_cm$CE_10_30)^2) / nrow(df_CE_10_30_cm))
rmse_CE_10_30 #quanto menor, melhor

rmse_CE_30_50 <- sqrt(sum((df_CE_30_50_cm$CE_30_50_reg - df_CE_30_50_cm$CE_30_50)^2) / nrow(df_CE_30_50_cm))
rmse_CE_30_50 #quanto menor, melhor



# preparado para krigar residuos -------------------------------------------------------
crs(CElab_train)<-crs

a_CE_0_10_reg = raster::extract(CE_0_10_pred, CElab_train, method="simple", sp=T); CE_0_10_reg = as.data.frame(a_CE_0_10_reg$CE_0_10_reg)

a_CE_10_30_reg = raster::extract(CE_10_30_pred, CElab_train, method="simple", sp=T);CE_10_30_reg = as.data.frame(a_CE_10_30_reg$CE_10_30_reg)

a_CE_30_50_reg = raster::extract(CE_30_50_pred, CElab_train, method="simple", sp=T); CE_30_50_reg = as.data.frame(a_CE_30_50_reg$CE_30_50_reg)

names(CE_0_10_reg)<-"CE_0_10_reg"
names(CE_10_30_reg)<-"CE_10_30_reg"
names(CE_30_50_reg)<-"CE_30_50_reg"


CE_reg<-cbind.data.frame(CElab_train,CE_0_10_reg,
                         CE_10_30_reg,CE_30_50_reg)
CE_reg

# CE_reg[, c(3:13)] <- sapply(CE_reg[, c(2:13)], as.numeric)

CE_reg$CE_0_10_res<-CE_reg$CE_0_10 - CE_reg$CE_0_10_reg
CE_reg$CE_10_30_res<-CE_reg$CE_10_30 - CE_reg$CE_10_30_reg
CE_reg$CE_30_50_res<-CE_reg$CE_30_50 - CE_reg$CE_30_50_reg

# krigando residuos -------------------------------------------------------

#0_10 cm
# id<-boxplot(CE_reg$CE_0_10_res,coef=2)
# id$stats
# lh <- quantile(CE_reg$CE_0_10_res,probs=0.25)#	Lower hinge (first quartile)
# uh <- quantile(CE_reg$CE_0_10_res,probs=0.75)	#Upper hinge (third quartile)
# step<- 1.5 * (uh-lh)	#Define the step as 1.5×IQR
# CE_reg$CE_0_10_res[CE_reg$CE_0_10_res < lh-step | CE_reg$CE_0_10_res > lh+step] <- NA
# CE_reg<-na.omit(CE_reg)

coordinates(CE_reg)<- ~X+Y
crs(CE_reg)<-crs

#0-10cm
variog.emp.CE_0_10_pred_res<- variogram(CE_0_10_res~1, CE_reg,cutoff = 300,boundaries=c(2:3*25, 2:11*50))
plot(variog.emp.CE_0_10_pred_res)
variog.prelim.CE_0_10_pred_res<- vgm(psill=600, model="Sph", range=300, nugget=200,boundaries=c(2:3*25, 2:11*50))
variog.aju.CE_0_10_pred_res <- fit.variogram(object=variog.emp.CE_0_10_pred_res, model=variog.prelim.CE_0_10_pred_res, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.CE_0_10_pred_res, variog.aju.CE_0_10_pred_res,main="CEa pred 0-10 cm")

krige.CE_0_10_pred_res<- krige(CE_0_10_res~1,locations=CE_reg, newdata=area_5, model=variog.aju.CE_0_10_pred_res)
spplot(krige.CE_0_10_pred_res["var1.pred"], main="CEa pred 0-10 cm")

#10-30 cm

variog.emp.CE_10_30_pred_res<- variogram(CE_10_30_res~1, CE_reg,cutoff = 229,boundaries=c(2:3*25, 2:11*50))
plot(variog.emp.CE_10_30_pred_res)
variog.prelim.CE_10_30_pred_res<- vgm(psill=60, model="Sph", range=150, nugget=35,boundaries=c(2:3*25, 2:11*50))
variog.aju.CE_10_30_pred_res <- fit.variogram(object=variog.emp.CE_10_30_pred_res, model=variog.prelim.CE_10_30_pred_res, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.CE_10_30_pred_res, variog.aju.CE_10_30_pred_res,main="CEa pred 10-30 cm")

krige.CE_10_30_pred_res<- krige(CE_10_30_res~1,locations=CE_reg, newdata=area_5, model=variog.aju.CE_10_30_pred_res)
spplot(krige.CE_10_30_pred_res["var1.pred"], main="CEa pred 10-30 cm")



#30-50 cm

variog.emp.CE_30_50_pred_res<- variogram(CE_30_50_res~1, CE_reg,cutoff = 250)
plot(variog.emp.CE_30_50_pred_res)
variog.prelim.CE_30_50_pred_res<- vgm(psill=40, model="Sph", range=180, nugget=30,boundaries=c(2:3*25, 2:11*50))
variog.aju.CE_30_50_pred_res <- fit.variogram(object=variog.emp.CE_30_50_pred_res, model=variog.prelim.CE_30_50_pred_res, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.CE_30_50_pred_res, variog.aju.CE_30_50_pred_res,main="CEa pred 30-50 cm")

krige.CE_30_50_pred_res<- krige(CE_30_50_res~1,locations=CE_reg, newdata=area_5, model=variog.aju.CE_30_50_pred_res)
spplot(krige.CE_30_50_pred_res["var1.pred"], main="CEa pred 30-50 cm")

# Somando residuos --------------------------------------------------------
#0_10 cm
CE_0_10_RK<-CE_0_10_pred[[1]] - krige.CE_0_10_pred_res[[1]]
CE_0_10_RK[CE_0_10_RK$CE_0_10_reg <= 0]<-0
spplot(CE_0_10_RK)

#10_30 cm
CE_10_30_RK<-CE_10_30_pred[[1]] - krige.CE_10_30_pred_res[[1]]
CE_10_30_RK[CE_10_30_RK$CE_10_30_reg <= 0]<-0
spplot(CE_10_30_RK)

#30_50 cm
CE_30_50_RK<-CE_30_50_pred[[1]] - krige.CE_30_50_pred_res[[1]]
CE_30_50_RK[CE_30_50_RK$CE_30_50_reg <= 0]<-0
spplot(CE_30_50_RK)

CE_pred_rk<-addLayer(CE_0_10_RK,
                     CE_10_30_RK,
                     CE_30_50_RK)

spplot(CE_pred_rk,col.regions=hcl.colors(20, palette = "Earth",rev=T),
       names.attr=c('(pred) EClab 0-10 cm','(pred) EClab10-30 cm','(pred) EClab 30-50 cm'))
