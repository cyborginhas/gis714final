#Load packages ----
require(biomod2) 
require(abind)
require(ade4)
require(caret)
require(checkmate)
require(dismo)
require(doParallel)
require(dplyr)
require(earth)
require(ecospat)
require(ENMeval)
require(foreach)
require(foreign)
require(gam)
require(gbm)
require(ggplot2)
require(Hmisc)
require(lattice)
require(MASS)
require(maxnet)
require(mda)
require(mgcv)
require(methods)
require(nnet)
require(parallel)
require(PresenceAbsence)
require(pROC)
require(purrr)
require(randomForest)
require(raster)
require(rasterVis)
require(reshape)
require(rlang)
require(rpart)
require(sp)
require(stats)
require(testthat)
require(tidyr)
require(utils)
require(rgdal)
setwd("~/Desktop/tohexports/")
#----
# 1. Setup study extent----
usa <- readOGR('~/Google Drive/Shared drives/Data/Vector/USA/us_lower_48_states.gpkg')
usa <- spTransform(usa, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
borders <- usa[usa$STATE_NAME%in%c('District of Columbia', 'Delaware', 'New Jersey', 'Maryland', 
                                   'West Virginia', 'Ohio','Pennsylvania','Virginia', 'New York', 'Kentucky','Tennessee'),]
# 2. Load species data----
aa.fia.ib <- read.csv2('~/Google Drive/My Drive/PhD/pops/slf/toh_exports/fia_ib.csv')[, c(4,5)] 
aa.fia.oob <- read.csv2('~/Google Drive/My Drive/PhD/pops/slf/toh_exports/fia_oob.csv')[, c(4,5)] 
names(aa.fia.ib) <- c('Latitude', 'Longitude')
aa.fia.ib<-aa.fia.ib[,c(2,1)]
aa.fia.ib.data<-aa.fia.ib
aa.fia.ib.data[,1:2]<-1
aa.fia.ib.pts <- SpatialPointsDataFrame(coords = aa.fia.ib, data=aa.fia.ib.data)
#plot(aa.fia.pts)
#plot(borders,add=TRUE)
crs(aa.fia.ib.pts) <-crs(borders)
aa.fia.ib.pts <- crop(aa.fia.ib.pts, borders)

# load the environmental raster layers (could be any supported format by the raster package)----
# Environmental variables extracted from Worldclim 
#myExpl_250m <- raster::getData('worldclim', download=T, var='bio', res=10)

biodir <- '~/Google Drive/Shared drives/APHIS  Projects/shared resources/data/worldclim1k/US/'
biovars <- stack(lapply(X=list.files(biodir), FUN=function(X){raster(paste(biodir, X, sep=''))}))
biovars <- stack(biovars[[1]],biovars[[12]])
e<-extent(borders)
e<-e*1.1
#roads.d <- raster('~/Google Drive/Shared drives/Data/Raster/Regional/roadsD_NE_1k.tif')
#pop<-raster('~/Google Drive/Shared drives/Data/Raster/Regional/popden_NE_1k.tif')
rails<-raster('~/Google Drive/Shared drives/Data/Raster/Regional/railsD_NE_1k.tif')
#canopy<-raster('~/Google Drive/Shared drives/Data/Raster/Regional/canop_NE_1k.tif')
pop<-raster('~/Google Drive/Shared drives/Data/Raster/Regional/popden_NE_250m.tif')
biovars<-crop(biovars,e)
rails<-crop(rails,e)
biovars <- resample(biovars, pop, method='bilinear')
rails <- resample(rails, pop, method='bilinear')

myExpl_250m <- mask(stack(biovars,rails), borders)
myExpl_250m <- crop(myExpl_250m, borders)
myExpl_250m <- stack(myExpl_250m)
names(myExpl_250m)<-c("BIO1","BIO12","rails")

#3. Rasterize response data to create presence and pseudoabsences----
#3A. FIA data
aa.fia.ras <- rasterize(x=aa.fia.pts, y=myExpl_250m[[1]], fun='count', background=0);
aa.fia.ras <- (aa.fia.ras*((myExpl_250m[[1]]*0)+1))>0
a2.fia.pts <- rasterToPoints(aa.fia.ras)
myRespName <- 'A_altissima_fia'
myResp <- a2.fia.pts[, 3] # the presence/absences data for our species
myResp[myResp==0] <- NA # setting 'true absences' to undefined
myRespXY <- a2.fia.pts[, c(1,2)] # the XY coordinates of species data
sum(na.omit(myResp))
myBiomodData.fia250m <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = myExpl_250m,
                                         resp.xy = myRespXY,
                                         resp.name = myRespName,
                                         PA.nb.rep = 1,
                                         PA.strategy = 'random',
                                         PA.nb.absences = sum(myResp, na.rm=T))
plot(myBiomodData.fia250m)
saveRDS(myBiomodData.fia250m, file = "myBiomodData.fia.rds")#144 presences

#3Ai. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()
#3Aii. Computing the models
getwd()#outputs placed here
myBiomodModelOut.fia250m <- biomod2::BIOMOD_Modeling(myBiomodData.fia250m,
                                                 models = c(#'CTA', 'SRE',  'MAXENT.Phillips','MAXENT.Phillips.2'
                                                   'GLM',
                                                   'GAM',
                                                   'MARS',
                                                   'FDA',
                                                   'GBM',
                                                   'RF',
                                                   'ANN'),
                                                 models.options = myBiomodOption,
                                                 NbRunEval = 5,
                                                 DataSplit = 100,
                                                 Prevalence = 0.5,
                                                 VarImport = 3,
                                                 models.eval.meth = 'TSS',
                                                 SaveObj = TRUE,
                                                 rescal.all.models = FALSE,
                                                 do.full.models = FALSE,
                                                 modeling.id=paste(myRespName,"FirstModeling.fia250m",sep=""))
saveRDS(myBiomodModelOut.fia250m, file = "myBiomodModelOut.fia250m.rds")
myBiomodModelEval.fia250m <- get_evaluations(myBiomodModelOut.fia250m) # get all models evaluation
dimnames(myBiomodModelEval.fia250m) # print the dimnames of this object
allmodels.fia250m.eval<-myBiomodModelEval.fia250m[c('TSS'),"Testing.data",,,] # print the eval scores of all selected models
allmodels.fia250m.eval<-as.data.frame(allmodels.fia250m.eval)
names(allmodels.fia250m.eval)<-"TSS"
models<-rownames(allmodels.fia250m.eval)
res<-250
db<-"fia"
allmodels<-cbind(models,as.data.frame(allmodels.fia250m.eval),res,db)
write.table(allmodels,"allmodels.csv", sep=",", row.names = FALSE,col.names = TRUE,append=TRUE)

vars_importance.fia250m <- as.data.frame(get_variables_importance(myBiomodModelOut.fia250m)) # print variable 
var_means<-rowMeans((vars_importance.fia250m),na.rm=TRUE)
res<-250
db<-"fia"
pred<-names(var_means)
var_means<-cbind(pred,as.data.frame(var_means),res,db)

#tiff('~/Desktop/tohexports/fia250mmeanvarimp.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
#barplot(colMeans(vars_imp_means_fia250m), ylab="Mean Variance Importance Score", xlab="Predictors")
#dev.off()
var_rank<-rank((vars_importance.fia250m)*-1,ties.method="average",na.last=NA)
var_rank<-cbind(pred,as.data.frame(var_rank),res,db)

write.table(var_means,"var_means.csv", row.names = FALSE, col.names = TRUE, append=TRUE, sep=",")#remove roads,canopy, and population
write.table(var_rank,"var_ranks.csv", row.names = FALSE, col.names = TRUE, append=TRUE, sep=",")#remove roads,canopy, and population

#3AiiiEnsembling the models
myBiomodEM.fia250m <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut.fia250m,
                                      chosen.models = 'all',
                                      em.by='all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.01),
                                      prob.mean = T,
                                      prob.cv = F, #don't use
                                      prob.ci = F, #prob.ci.alpha = 0.05,
                                      prob.median = F,
                                      committee.averaging = F,
                                      prob.mean.weight = T,
                                      prob.mean.weight.decay = 'proportional' )

saveRDS(myBiomodEM.fia250m, file = "myBiomodEM.fia250m.rds")
ensmodeleval<-get_evaluations(myBiomodEM.fia250m) # get evaluation scores
ensmodeleval<-ensmodeleval[[2]]
stats<-row.names(ensmodeleval)
ensmodeleval<-cbind(stats,as.data.frame(ensmodeleval))
write.table(append = T, sep=",", ensmodeleval, "enseval.csv", row.names = FALSE, col.names = TRUE)

###3Aiv. projection over the globe under current conditions
myBiomodProj.fia250m <- BIOMOD_Projection(modeling.output = myBiomodModelOut.fia250m,
                                  new.env = myExpl_250m,
                                  proj.name = 'current.fia',
                                  selected.models = 'all',
                                  binary.meth = 'TSS',
                                  compress = 'xz',
                                  clamping.mask = F,
                                  output.format = '.grd')
saveRDS(myBiomodProj.fia250m, file = "myBiomodProj.fia250m.rds")
myCurrentProj.fia <- get_predictions(myBiomodProj.fia250m) # if you want to make custom plots, you can also get the projected map
myBiomodEF.fia250m <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM.fia250m,
                                         projection.output = myBiomodProj.fia250m)
saveRDS(myBiomodEF.fia250m, file = "myBiomodEF.fia250m.rds")
plot(myBiomodEF.fia250m)
pred.out.fia250m <- myBiomodEF.fia250m@proj@val[[1]]
pred.out.fia250m[]<-pred.out.fia250m[]/1000
plot(pred.out.fia250m)
writeRaster(pred.out.fia250m, filename = 'toh_ens_fia250m_wgs.tif', format="GTiff")
crs<-"+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" 
pred.out.fia250m.lcc<-(projectRaster(pred.out.fia250m,crs=crs))
writeRaster(pred.out.fia250m.lcc, filename = 'toh_ens_fia250m_lcc.tif', format="GTiff")

#3Av. output plot
rpts <- rasterToPoints(pred.out.fia250m.lcc)
rdf <- as.data.frame(rpts)
ggsdm <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
  geom_path(data=borders, aes(x=long, y=lat, group=group), col='white', lwd=1.1, alpha=.3) + 
  scale_fill_continuous(type='viridis') + theme_void() + theme(legend.position='none')
png(paste('pred.out.fia.250m.lcc', 
          gsub(':', '', substr(Sys.time(), 12, 19)), '.png', sep=''), 
    height=1080, width=2160); plot(ggsdm); dev.off()

