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
# load("D:/OneDrive/Revistas/Soil Systems/salinidade/data/mapas_CE_SM_DIBA_log.RData")

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

plot(CElab_train)
points(CElab_valid,col="red")

# mapa de bolhas ----------------------------------------------------------

bub_CE_0_10<-bubble(CElab["CE_0_10"])
bub_CE_10_30<-bubble(CElab["CE_10_30"])
bub_CE_30_50<-bubble(CElab["CE_30_50"])

bub_CE_lab<-c(bub_CE_0_10,
              bub_CE_10_30,
              bub_CE_30_50,x.same = NA, y.same = NA,
              layout = NULL, merge.legends = FALSE, recursive = FALSE)

bub_CE_lab



# Krig EM38 nos dados de CE lab ---------------------------------------------------------------
coordinates(EM38)<-~X_UTM_24S+Y_UTM_24S
EM38<-EM38[-zerodist(EM38)[,1],]
crs(EM38)<-crs

set.seed(1)
EM_38_size<- floor(0.95 * nrow(EM38)) #25%
set.seed(123)
EM38_retirar <- sample(seq_len(nrow(EM38)), size = EM_38_size)
EM38_2 <- EM38[-EM38_retirar, ] # conjunto de dados - 25%

tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/hist_EM38.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

graphics::hist(as.data.frame(EM38_2)[9:12])

dev.off()

tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/hist_EM38_log.tiff", height = 15, width = 15, units="cm",compression = "lzw", res = 300)

graphics::hist(log(as.data.frame(EM38_2)[9:12]))

dev.off()

crs(EM38_2)<-crs
crs(area_5)<-crs


par(mfrow=(c(1,2)))


plot(EM38,main="100% ou 20669 pontos",pch=10,cex=0.2)
plot(EM38_2, main= "5% ou 1034 pontos",pch=10,cex=0.2)

# CE 1 m
variog.emp.EM_1m<- variogram(log(CV_1_0m) ~ 1, EM38_2, cutoff = 90)
plot(variog.emp.EM_1m)
variog.prelim.EM_1m<- vgm(psill=0.2, model="Sph", range=85, nugget=0.02,boundaries=c(2:3*25, 2:11*50))
variog.aju.EM_1m <- fit.variogram(object=variog.emp.EM_1m, model=variog.prelim.EM_1m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.EM_1m, variog.aju.EM_1m,main="CEa bobina de 1 m")

krige.CE_1m_lab<- krigeTg(CV_1_0m~1,locations=EM38_2, newdata=CElab, model=variog.aju.EM_1m,lambda=0)
spplot(krige.CE_1m_lab["var1TG.pred"], main="CEa bobina 1 metro ")

#CE0.5
variog.emp.EM_0_5m<- variogram(log(CV_0_5m) ~ 1, EM38_2, cutoff = 90)
plot(variog.emp.EM_0_5m)
variog.prelim.EM_0_5m<- vgm(psill=0.42, model="Sph", range=85, nugget=0.05)
variog.aju.EM_0_5m <- fit.variogram(object=variog.emp.EM_0_5m, model=variog.prelim.EM_0_5m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.EM_0_5m, variog.aju.EM_0_5m,main="CEa bobina de 0.5 m")

krige.CE_0_5m_lab<- krigeTg(CV_0_5m~1,locations=EM38_2, newdata=CElab, model=variog.aju.EM_0_5m,lambda=0)
spplot(krige.CE_0_5m_lab["var1TG.pred"], main="CEa bobina 0.5 metro ")

#SM1m
variog.emp.SM_1m<- variogram(log(IV_1_0m) ~ 1, EM38_2, cutoff = 180)
plot(variog.emp.SM_1m)
variog.prelim.SM_1m<- vgm(psill=0.15, model="Sph", range=160, nugget=0.005)
variog.aju.SM_1m <- fit.variogram(object=variog.emp.SM_1m, model=variog.prelim.SM_1m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.SM_1m, variog.aju.SM_1m,main="SMa bobina de 1m")

krige.SM_1m_lab<- krigeTg(IV_1_0m~1,locations=EM38_2, newdata=CElab, model=variog.aju.SM_1m,lambda=0)
spplot(krige.SM_1m_lab["var1TG.pred"], main="SMa bobina 1 metro ")

#SM0.5
variog.emp.SM_0_5m<- variogram(log(IV_0_5m) ~ 1, EM38_2, cutoff = 180)
plot(variog.emp.SM_0_5m)
variog.prelim.SM_0_5m<- vgm(psill=0.9, model="Sph", range=150, nugget=0.28)
variog.aju.SM_0_5m <- fit.variogram(object=variog.emp.SM_0_5m, model=variog.prelim.SM_0_5m, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.SM_0_5m, variog.aju.SM_0_5m,main="SMa bobina de 0.5 m")

krige.SM_0_5m_lab<- krigeTg(IV_0_5m~1,locations=EM38_2, newdata=CElab, model=variog.aju.SM_0_5m,lambda=0)
spplot(krige.SM_0_5m_lab["var1TG.pred"], main="SMa bobina 0.5 metro ")

# plotando semivariogramas juntos -----------------------------------------

par(mfrow=c(2,2))

plot(variog.emp.EM_1m$dist, variog.emp.EM_1m$gamma, xlab="Distance",
     ylab="Semivariance",main="aEC (coil 1 m)")

lines(variogramLine(variog.aju.EM_1m, maxdist=100))


plot(variog.emp.EM_0_5m$dist, variog.emp.EM_0_5m$gamma, xlab="Distance",
     ylab="Semivariance",main="aEC (coil 0.5 m)")

lines(variogramLine(variog.aju.EM_0_5m, maxdist=100))

plot(variog.emp.SM_1m$dist, variog.emp.SM_1m$gamma, xlab="Distance",
     ylab="Semivariance",main="aMS (coil 1 m)")

lines(variogramLine(variog.aju.SM_1m, maxdist=200))

plot(variog.emp.SM_0_5m$dist, variog.emp.SM_0_5m$gamma, xlab="Distance",
     ylab="Semivariance",main="aMS (coil 0.5 m)")

lines(variogramLine(variog.aju.SM_0_5m, maxdist=200))


# Convertendo dados krigado em dataframe ----------------------------------
CE_1m_TG_lab<-data.frame(krige.CE_1m_lab["var1TG.pred"])
names(CE_1m_TG_lab)<-c("CV_1_0m","X","Y","opt")

CE_0_5_m_TG_lab<-data.frame(krige.CE_0_5m_lab["var1TG.pred"])
names(CE_0_5_m_TG_lab)<-c("CV_0_5m","X","Y","opt")

SM_1_m_TG_lab<-data.frame(krige.SM_1m_lab["var1TG.pred"])
names(SM_1_m_TG_lab)<-c("IV_1_0m","X","Y","opt")

SM_0_5_m_TG_lab<-data.frame(krige.SM_0_5m_lab["var1TG.pred"])
names(SM_0_5_m_TG_lab)<-c("IV_0_5m","X","Y","opt")

sensor_to_CE_points<-bind_cols(as.data.frame(CElab),CE_1m_TG_lab[1], CE_0_5_m_TG_lab[1],
                               SM_1_m_TG_lab[1],SM_0_5_m_TG_lab[1])
sensor_to_CE_points

# Regressoes --------------------------------------------------------------
#CE 0-10
EC_0_10_lab<-regsubsets(CE_0_10 ~ .,data=sensor_to_CE_points[,c(4,10:13)],nbest = 50, nvmax = 50,method = "exhaustive",really.big=T)

EC_0_10_reg_lab<- summary(EC_0_10_lab)

coef(EC_0_10_lab, id=which(EC_0_10_reg_lab$bic == min(EC_0_10_reg_lab$bic)))
coef(EC_0_10_lab, id=which(EC_0_10_reg_lab$adjr2 == max(EC_0_10_reg_lab$adjr2)))

EC_0_10_lm_r2_lab<-lm(CE_0_10 ~  CV_0_5m,data=sensor_to_CE_points)
summary(EC_0_10_lm_r2_lab)

#CE 10-30
EC_10_30_lab<-regsubsets(CE_10_30 ~ .,data=sensor_to_CE_points[,c(5,10:13)],nbest = 50, nvmax = 50,method = "exhaustive",really.big=T)

EC_10_30_reg_lab<- summary(EC_10_30_lab)

coef(EC_10_30_lab, id=which(EC_10_30_reg_lab$bic == min(EC_10_30_reg_lab$bic)))
coef(EC_10_30_lab, id=which(EC_10_30_reg_lab$adjr2 == max(EC_10_30_reg_lab$adjr2)))

EC_10_30_lm_r2_lab<-lm(CE_10_30 ~ CV_1_0m+CV_0_5m+IV_1_0m+IV_0_5m ,data=sensor_to_CE_points)
summary(EC_10_30_lm_r2_lab)

#CE 10-30
EC_30_50_lab<-regsubsets(CE_30_50 ~ .,data=sensor_to_CE_points[,c(6,10:13)],nbest = 50, nvmax = 50,method = "exhaustive",really.big=T)
EC_30_50_reg_lab<- summary(EC_30_50_lab)

coef(EC_30_50_lab, id=which(EC_30_50_reg_lab$bic == min(EC_30_50_reg_lab$bic)))
coef(EC_30_50_lab, id=which(EC_30_50_reg_lab$adjr2 == max(EC_30_50_reg_lab$adjr2)))

EC_30_50_lm_r2_lab<-lm(CE_30_50 ~ CV_1_0m + CV_0_5m+IV_1_0m ,data=sensor_to_CE_points)
summary(EC_30_50_lm_r2_lab)


tab_model(
        EC_0_10_lm_r2_lab,EC_10_30_lm_r2_lab,EC_30_50_lm_r2_lab,
        dv.labels = c("CE (pred 0-10 cm) ", "CE (pred 10-30 cm)","CE (pred 30-50 cm)"),
        string.pred = "Coeffcient",
        string.ci = "Conf. Int (95%)",
        string.p = "P-Value"
)


# predict -----------------------------------------------------------------
CE_0_10_pred_lab<-predict(EC_0_10_lm_r2_lab,EM38_2)
CE_10_30_pred_lab<-predict(EC_10_30_lm_r2_lab,EM38_2)
CE_30_50_pred_lab<-predict(EC_30_50_lm_r2_lab,EM38_2)


# analise exploratoria dos dados de CElab pred ----------------------------
# removendo todos os dados negativos de CE.
CE_0_10_pred_lab<-as.data.frame(CE_0_10_pred_lab)
CE_0_10_pred_lab$X<-data.frame(EM38_2@coords)[,1]
CE_0_10_pred_lab$Y<-data.frame(EM38_2@coords)[,2]

CE_10_30_pred_lab<-as.data.frame(CE_10_30_pred_lab)
CE_10_30_pred_lab$X<-data.frame(EM38_2@coords)[,1]
CE_10_30_pred_lab$Y<-data.frame(EM38_2@coords)[,2]

CE_30_50_pred_lab<-as.data.frame(CE_30_50_pred_lab)
CE_30_50_pred_lab$X<-data.frame(EM38_2@coords)[,1]
CE_30_50_pred_lab$Y<-data.frame(EM38_2@coords)[,2]


CE_0_10_pred_lab[CE_0_10_pred_lab <= 0] <- NA
CE_10_30_pred_lab[CE_10_30_pred_lab <= 0] <- NA
CE_30_50_pred_lab[CE_30_50_pred_lab <= 0] <- NA

CE_0_10_pred_lab<-na.omit(CE_0_10_pred_lab)
CE_10_30_pred_lab<-na.omit(CE_10_30_pred_lab)
CE_30_50_pred_lab<-na.omit(CE_30_50_pred_lab)


hist(CE_0_10_pred_lab$CE_0_10_pred_lab)
hist(CE_10_30_pred_lab$CE_10_30_pred_lab)
hist(CE_30_50_pred_lab$CE_30_50_pred_lab)

# krigando CElab pred por EM38 --------------------------------------------
coordinates(CE_0_10_pred_lab)<-~X+Y
coordinates(CE_10_30_pred_lab)<-~X+Y
coordinates(CE_30_50_pred_lab)<-~X+Y

crs(CE_0_10_pred_lab)<-crs
crs(CE_10_30_pred_lab)<-crs
crs(CE_30_50_pred_lab)<-crs


# CE 0_10 cm
variog.emp.CE_0_10_pred<- variogram(CE_0_10_pred_lab ~ 1, CE_0_10_pred_lab, cutoff = 80)
plot(variog.emp.CE_0_10_pred)
variog.prelim.CE_0_10_pred<- vgm(psill=26, model="Sph", range=75, nugget=2.9,boundaries=c(2:3*25, 2:11*50))
variog.aju.CE_0_10_pred <- fit.variogram(object=variog.emp.CE_0_10_pred, model=variog.prelim.CE_0_10_pred, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.CE_0_10_pred, variog.aju.CE_0_10_pred,main="CEa bobina de 1 m")

# krige.CE_0_10_pred<- krige(CE_0_10_pred_lab~1,locations=CE_0_10_pred_lab, newdata=area_5, model=variog.aju.CE_0_10_pred)
# spplot(krige.CE_0_10_pred["var1.pred"], main="CEa 0-10 cm predita por EM38")

raster::writeRaster(raster(krige.CE_0_10_pred["var1.pred"]), filename="D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_0_10_pred.tif", overwrite=TRUE,bylayer=T)
CE_0_10_pred_pontos<-raster("D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_0_10_pred.tif")


# CE 10_30 cm
variog.emp.CE_10_30_pred<- variogram(CE_10_30_pred_lab ~ 1, CE_10_30_pred_lab, cutoff = 180)
plot(variog.emp.CE_10_30_pred)
variog.prelim.CE_10_30_pred<- vgm(psill=38, model="Sph", range=140, nugget=14,boundaries=c(2:3*25, 2:11*50))
variog.aju.CE_10_30_pred <- fit.variogram(object=variog.emp.CE_10_30_pred, model=variog.prelim.CE_10_30_pred, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.CE_10_30_pred, variog.aju.CE_10_30_pred,main="CEa bobina de 1 m")

# krige.CE_10_30_pred<- krige(CE_10_30_pred_lab~1,locations=CE_10_30_pred_lab, newdata=area_5, model=variog.aju.CE_10_30_pred)
# spplot(krige.CE_10_30_pred["var1.pred"], main="CEa 10-30 cm predita por EM38")

# raster::writeRaster(raster(krige.CE_10_30_pred["var1.pred"]), filename="D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_10_30_pred.tif", overwrite=TRUE,bylayer=T)
CE_10_30_pred_pontos<-raster("D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_10_30_pred.tif")


# CE 30_50 cm
variog.emp.CE_30_50_pred<- variogram(CE_30_50_pred_lab ~ 1, CE_30_50_pred_lab, cutoff = 130)
plot(variog.emp.CE_30_50_pred)
variog.prelim.CE_30_50_pred<- vgm(psill=32, model="Sph", range=110, nugget=4,boundaries=c(2:3*25, 2:11*50))
variog.aju.CE_30_50_pred <- fit.variogram(object=variog.emp.CE_30_50_pred, model=variog.prelim.CE_30_50_pred, fit.sills=F, fit.ranges=F, fit.method = 7)
plot(variog.emp.CE_30_50_pred, variog.aju.CE_30_50_pred,main="CEa bobina de 1 m")

# krige.CE_30_50_pred<- krige(CE_30_50_pred_lab~1,locations=CE_30_50_pred_lab, newdata=area_5, model=variog.aju.CE_30_50_pred)
spplot(krige.CE_30_50_pred["var1.pred"], main="CEa 30-50 cm predita por EM38")

# raster::writeRaster(raster(krige.CE_30_50_pred["var1.pred"]), filename="D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_30_50_pred.tif", overwrite=TRUE,bylayer=T)
CE_30_50_pred_pontos<-raster("D:/OneDrive/Revistas/Soil Systems/salinidade/data/CE_30_50_pred.tif")



# plotando semivariogramas juntos krigagem CE lab predita -----------------------------------------

tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/semivariograma_CE_lab_pred_abordagem_2.tiff", height = 7, width = 18, units="cm",compression = "lzw", res = 300)

par(mfrow=c(1,3))

plot(variog.emp.CE_0_10_pred$dist, variog.emp.CE_0_10_pred$gamma, xlab="Distance",
     ylab="Semivariance",main="pred aEC 0_10 cm")

lines(variogramLine(variog.aju.CE_0_10_pred, maxdist=100))


plot(variog.emp.CE_10_30_pred$dist, variog.emp.CE_10_30_pred$gamma, xlab="Distance",
     ylab="Semivariance",main="pred aEC lab 10_30 cm")

lines(variogramLine(variog.aju.CE_10_30_pred, maxdist=180))

plot(variog.emp.CE_30_50_pred$dist, variog.emp.CE_30_50_pred$gamma, xlab="Distance",
     ylab="Semivariance",main="pred aEC lab 30_50 cm")

lines(variogramLine(variog.aju.CE_30_50_pred, maxdist=200))

dev.off()



# plotando mapas de CElab pred juntos -------------------------------------

crs(CE_0_10_pred_pontos)<-"+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
CE_0_10_pred_pontos[CE_0_10_pred_pontos$CE_0_10_pred <= 0]<-0
spplot(CE_0_10_pred_pontos)

crs(CE_10_30_pred_pontos)<-"+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
CE_10_30_pred_pontos[CE_10_30_pred_pontos$CE_10_30_pred <= 0]<-0
spplot(CE_10_30_pred_pontos)

crs(CE_30_50_pred_pontos)<-"+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
CE_30_50_pred_pontos[CE_30_50_pred_pontos$CE_30_50_pred <= 0]<-0
spplot(CE_30_50_pred_pontos)


CElab_pred_pontos<-addLayer(CE_0_10_pred_pontos,
                            CE_10_30_pred_pontos,
                            CE_30_50_pred_pontos)

tiff("D:/OneDrive/Revistas/Soil Systems/salinidade/results/images/EC_pred_abordagem2.tiff", height = 10, width = 18, units="cm",compression = "lzw", res = 300)

spplot(CElab_pred_pontos,col.regions=hcl.colors(50, palette = "Earth",rev=T),
       names.attr=c('(pred) EC lab 0-10 cm','(pred) EC lab10-30 cm','(pred) EC lab 30-50 cm'),
       at= seq(0, 45, 1))

dev.off()


# Preparando para validacao ---------------------------------------------------------------

a_CE_0_10_cm = raster::extract(CE_0_10_pred_pontos, CElab, method="simple", sp=T); df_CE_0_10_cm = as.data.frame(a_CE_0_10_cm)
a_CE_10_30_cm = raster::extract(CE_10_30_pred_pontos, CElab, method="simple", sp=T); df_CE_10_30_cm = as.data.frame(a_CE_10_30_cm)
a_CE_30_50_cm = raster::extract(CE_30_50_pred_pontos, CElab, method="simple", sp=T); df_CE_30_50_cm = as.data.frame(a_CE_30_50_cm)


# calculando erros --------------------------------------------------------

res_CE_0_10_ponto<-mean(df_CE_0_10_cm$CE_0_10 - df_CE_0_10_cm$CE_0_10_pred)
res_CE_0_10_ponto

res_CE_10_30_ponto<-mean(df_CE_10_30_cm$CE_10_30 - df_CE_10_30_cm$CE_10_30_pred)
res_CE_10_30_ponto

res_CE_30_50_ponto<-mean(df_CE_30_50_cm$CE_30_50 - df_CE_30_50_cm$CE_30_50_pred)
res_CE_30_50_ponto

rmse_CE_0_10_ponto <- sqrt(sum((df_CE_0_10_cm$CE_0_10_pred - df_CE_0_10_cm$CE_0_10)^2) / nrow(df_CE_0_10_cm))
rmse_CE_0_10_ponto #quanto menor, melhor

rmse_CE_10_30_ponto <- sqrt(sum((df_CE_10_30_cm$CE_10_30_pred - df_CE_10_30_cm$CE_10_30)^2) / nrow(df_CE_10_30_cm))
rmse_CE_10_30_ponto #quanto menor, melhor

rmse_CE_30_50_ponto <- sqrt(sum((df_CE_30_50_cm$CE_30_50_pred - df_CE_30_50_cm$CE_30_50)^2) / nrow(df_CE_30_50_cm))
rmse_CE_30_50_ponto #quanto menor, melhor

