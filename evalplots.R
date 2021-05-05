#load libraries
library(data.table)#clean, manipulate data
library(ggplot2)#plotting
library(raster)
library(rgdal)
library(ecospat)

#ensemble TSS
enseval<-read.csv("~/Desktop/tohexports/enseval.csv")
enseval<-as.data.table(enseval)
TSSens<-enseval[stats=="TSS",]
TSSens$res=c(1000,250,1000,250)
TSSens$db=c("fia","cs","cs","fia")
TSS<-TSSens[,.(models="*ENS",TSS=as.numeric(Testing.data),res,db)]
TSSens_Sens<-TSSens[,.(Percent=round(as.numeric(Sensitivity),1),res,db,Type="sensitivity")]
TSSens_Spec<-TSSens[,.(Percent=round(as.numeric(Specificity),1),res,db,Type="specificity")]
TSSens_t<-rbind(TSSens_Sens,TSSens_Spec)

#TSS all models
alleval<-read.csv("~/Desktop/tohexports/allmodels.csv")
alleval<-as.data.table(alleval)
alleval<-alleval[models!="models",.(models,TSS=as.numeric(TSS), res=as.numeric(res),db)]
allens_eval<-rbind(alleval,TSS)
allens_eval<-as.data.frame(allens_eval)
allens_eval$db<-as.factor(allens_eval$db)
allens_eval$models<-as.factor(allens_eval$models)
levels(allens_eval$db)<-c("haphazard","systematic")

#Plot TSS for all models & export plot as png
p<-ggplot(allens_eval, aes(fill=as.factor(res), y=TSS, x=models)) + 
  geom_bar(position="dodge", stat="identity", lwd=0) +
  scale_fill_manual(values=c("grey50","grey85"), name="Grain size (m)") +
  facet_wrap(~db) + theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.position = c(.15,.88), legend.direction="vertical", 
        legend.background = element_blank(),axis.text.x = element_text(angle = 90)) + labs(x="")

png("~/Desktop/allmodels.png", units="in", width=8, height=6, res=300)
p
dev.off()

#Plot TSS for ensemble models & export plot as png
TSSens_t<-as.data.frame(TSSens_t)
TSSens_t$db<-as.factor(TSSens_t$db)
levels(TSSens_t$db)<-c("haphazard","systematic")

p2<-ggplot(TSSens_t, aes(fill=as.factor(res), y=Percent, x=Type)) + 
  geom_bar(position="dodge", stat="identity", lwd=0) +
  scale_fill_manual(values=c("grey50","grey85"), name="Grain size (m)") +
  facet_wrap(~db) + theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        legend.position="none") + labs(x="") + ylim(0,100) + labs(y="Mean Performance Value (%)") + 
  geom_text(aes(label=Percent, y=Percent+3), position = position_dodge(0.9))
png("~/Desktop/ensmodels.png", units="in", width=8, height=6, res=300)
p2
dev.off()

#Variable importance - data manipulation
var_ranks<-read.csv("~/Desktop/tohexports/var_ranks.csv")
var_ranks<-as.data.table(var_ranks)
var_ranks<-var_ranks[pred!="pred",.(pred=pred,rank=as.numeric(var_rank), res=as.factor(res), db=as.factor(db))]
var_ranks[,.(median(rank)), by=.(pred)]
var_means<-read.csv("~/Desktop/tohexports/var_means.csv")
var_means<-as.data.table(var_means)
var_means<-var_means[pred!="pred",.(pred=pred,mean=as.numeric(var_means), res=as.factor(res), db=as.factor(db))]
var_means[,.(mean(mean)), by=.(pred)]

#Testing data evaluation
fia_testing<-readOGR("~/Desktop/tohoutputs/aa.fia.oob.lcc.shp")
fia_1k<-raster("~/Desktop/tohoutputs/toh_ens_fia1k_lcc.tif")
fia_1k_boyce<-ecospat.boyce(fia_1k, coordinates(fia_testing), nclass=0, window.w="default", res=100, PEplot = TRUE)
fia_1k_boyce_t<-as.data.frame(cbind(fia_1k_boyce$HS,fia_1k_boyce$F.ratio))
names(fia_1k_boyce_t)<-c("Suitability","Pred/Exp (F) Ratio")

png("~/Desktop/fia1kboyce.png", units="in", width=8, height=6, res=300)
plot(fia_1k_boyce_t, col="white", cex.lab=1.75, cex.axis=1.5,mgp=c(2.5,0.75,0))
points(fia_1k_boyce_t[fia_1k_boyce_t$`Pred/Exp (F) Ratio`>=1,], col="#01665e",pch=20,lwd=3)
points(fia_1k_boyce_t[fia_1k_boyce_t$`Pred/Exp (F) Ratio`<=1,], col="#80cdc1",pch=20)
text.default(0.63, 3, labels="threshold = 0.64", srt=90,offset=0.5, cex=1.25)
abline(v=0.6424195)
dev.off()


#fia_250<-raster("~/Desktop/tohoutputs/toh_ens_fia250m_lcc.tif")
fia_250_boyce<-ecospat.boyce(fia_250, coordinates(fia_testing), nclass=0, window.w="default", res=100, PEplot = TRUE)
fia_250_boyce_t<-as.data.frame(cbind(fia_250_boyce$HS,fia_250_boyce$F.ratio))
names(fia_250_boyce_t)<-c("Suitability","Pred/Exp (F) Ratio")
png("~/Desktop/fia250boyce.png", units="in", width=8, height=6, res=300)
plot(fia_250_boyce_t, col="white", cex.lab=1.75, cex.axis=1.5,mgp=c(2.5,0.75,0))
points(fia_250_boyce_t[fia_250_boyce_t$`Pred/Exp (F) Ratio`>=1,], col="#01665e",pch=20,lwd=3)
points(fia_250_boyce_t[fia_250_boyce_t$`Pred/Exp (F) Ratio`<=1,], col="#80cdc1",pch=20)
text.default(0.52, 3.25, labels="threshold = 0.53", srt=90,offset=0.5, cex=1.25)
abline(v=0.5329216)
dev.off()

cs_testing<-readOGR("~/Desktop/tohoutputs/aa.cs.oob.lcc.shp")
cs_1k<-raster("~/Desktop/tohoutputs/toh_ens_cs1k_lcc.tif")
cs_1k_boyce<-ecospat.boyce(cs_1k, coordinates(cs_testing), nclass=0, window.w="default", res=100, PEplot = TRUE)
cs_1k_boyce_t<-as.data.frame(cbind(cs_1k_boyce$HS,cs_1k_boyce$F.ratio))
names(cs_1k_boyce_t)<-c("Suitability","Pred/Exp (F) Ratio")
png("~/Desktop/cs1kboyce.png", units="in", width=8, height=6, res=300)
plot(cs_1k_boyce_t, col="white", cex.lab=1.75, cex.axis=1.5,mgp=c(2.5,0.75,0))
points(cs_1k_boyce_t[cs_1k_boyce_t$`Pred/Exp (F) Ratio`>=1,], col="#01665e",pch=20,lwd=3)
points(cs_1k_boyce_t[cs_1k_boyce_t$`Pred/Exp (F) Ratio`<=1,], col="#80cdc1",pch=20)
text.default(0.65, 5, labels="threshold = 0.66", srt=90,offset=0.5, cex=1.25)
abline(v=0.66417292)
dev.off()

cs_250m<-raster("~/Desktop/tohoutputs/toh_ens_cs250m_lcc.tif")
#cs_250m_boyce<-ecospat.boyce(cs_250m, coordinates(cs_testing), nclass=0, window.w="default", res=100, PEplot = TRUE)
cs_250m_boyce_t<-as.data.frame(cbind(cs_250m_boyce$HS,cs_250m_boyce$F.ratio))
names(cs_250m_boyce_t)<-c("Suitability","Pred/Exp (F) Ratio")
png("~/Desktop/cs250mboyce.png", units="in", width=8, height=6, res=300)
plot(cs_250m_boyce_t, col="white", cex.lab=1.75, cex.axis=1.5,mgp=c(2.5,0.75,0))
points(cs_250m_boyce_t[cs_250m_boyce_t$`Pred/Exp (F) Ratio`>=1,], col="#01665e",pch=20,lwd=3)
points(cs_250m_boyce_t[cs_250m_boyce_t$`Pred/Exp (F) Ratio`<=1,], col="#80cdc1",pch=20)
text.default(0.62, 7, labels="threshold = 0.63", srt=90,offset=0.5, cex=1.25)
abline(v=0.63427757)
dev.off()

fia1k_bn<-fia_1k
fia1k_bn[fia1k_bn[]>=0.6424195]<-1
fia1k_bn[fia1k_bn[]<0.6424195]<-0
plot(fia1k_bn)
fia1kcorr<-raster::extract(fia1k_bn, fia_testing, method="simple")
fia1kcorr<-na.omit(sum(fia1kcorr)/length(fia1kcorr))
fia1kcorr_wcs<-raster::extract(fia1k_bn, cs_testing, method="simple")
fia1kcorr_wcs<-fia1kcorr_wcs[!is.na(fia1kcorr_wcs)]
fia1kcorr_wcs<-sum(fia1kcorr_wcs)/length(fia1kcorr_wcs)


fia250_bn<-fia_250
fia250_bn[fia250_bn[]>=0.5329216]<-1
fia250_bn[fia250_bn[]<0.5329216]<-0
plot(fia250_bn)
fia250corr<-raster::extract(fia250_bn, fia_testing, method="simple")
fia250corr<-na.omit(sum(fia250corr)/length(fia250corr))
fia250corr_wcs<-raster::extract(fia250_bn, cs_testing, method="simple")
fia250corr_wcs<-fia250corr_wcs[!is.na(fia250corr_wcs)]
fia250corr_wcs<-sum(fia250corr_wcs)/length(fia250corr_wcs)


cs1k_bn<-cs_1k
cs1k_bn[cs1k_bn[]>=0.66417292]<-1
cs1k_bn[cs1k_bn[]<0.66417292]<-0
plot(cs1k_bn)
cs1kcorr<-raster::extract(cs1k_bn, cs_testing, method="simple")
cs1kcorr<-cs1kcorr[!is.na(cs1kcorr)]
cs1kcorr<-(sum(cs1kcorr)/length(cs1kcorr))
cs1kcorr_wfia<-raster::extract(cs1k_bn, fia_testing, method="simple")
cs1kcorr_wfia<-cs1kcorr_wfia[!is.na(cs1kcorr_wfia)]
cs1kcorr_wfia<-(sum(cs1kcorr_wfia)/length(cs1kcorr_wfia))


cs250_bn<-cs_250m
cs250_bn[cs250_bn[]>=0.63427757]<-1
cs250_bn[cs250_bn[]<0.63427757]<-0
plot(cs250_bn)
cs250corr<-raster::extract(cs250_bn, cs_testing, method="simple")
cs250corr<-cs250corr[!is.na(cs250corr)]
cs250corr<-(sum(cs250corr)/length(cs250corr))
cs250corr_wfia<-raster::extract(cs250_bn, fia_testing, method="simple")
cs250corr_wfia<-cs250corr_wfia[!is.na(cs250corr_wfia)]
cs250corr_wfia<-(sum(cs250corr_wfia)/length(cs250corr_wfia))

writeRaster(fia1k_bn,"~/Desktop/tohoutputs/fia1kbn.tif", format="GTiff")
writeRaster(fia250_bn,"~/Desktop/tohoutputs/fia250bn.tif", format="GTiff")
writeRaster(cs1k_bn,"~/Desktop/tohoutputs/cs1kbn.tif", format="GTiff")
writeRaster(cs250_bn,"~/Desktop/tohoutputs/cs250bn.tif", format="GTiff")
