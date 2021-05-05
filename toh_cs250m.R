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
#FIA data----
#FIA pts plot
#tiff('fiapts.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
#plot(aa.fia.pts,pch=21, axes=TRUE,bg='blue', xlim=c(-90,-72))
#plot(borders, add=TRUE)
#dev.off()

#Citizen science data----
aa.cs.ib <- read.csv2('~/Google Drive/My Drive/PhD/pops/slf/toh_exports/cs_ib.csv')[, c(4,5)]
aa.cs.oob <- read.csv2('~/Google Drive/My Drive/PhD/pops/slf/toh_exports/cs_oob.csv')[, c(4,5)] 
names(aa.cs.ib) <- c('Latitude', 'Longitude')
aa.cs.ib<-aa.cs.ib[,c(2,1)]
aa.cs.pts <- SpatialPoints(coords = aa.cs.ib)
crs(aa.cs.pts) <-crs(borders)
aa.cs.pts <- crop(aa.cs.pts, borders)
plot(aa.cs.pts)
#CS pts plot
#tiff('fiapts.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
#plot(aa.cs.pts,pch=21, axes=TRUE,bg='salmon', xlim=c(-90,-72))
#plot(borders, add=TRUE)
#dev.off()

#Examine distribution of species data points (FIA vs. CS)----
#Create hex grids of various sizes across study area, pull out areas where fia sampling effort same as cs
#cs_l<-seq(from = 0.2, to = 2, by=0.1)
#list_id<-seq(from=1, to= 19, by=1)
#samps<-as.data.frame(cbind(cs_l,list_id))
#my_list<-list()#list of sample point data derived from individual hex polygons (by res)
#my_list2<-list()#list of hex polygons associated individual hex polygons (by res) 
#library(GISTools)
#borders<-gUnion(borders,borders)
#plot(borders)
#for(i in 1:length(samps$cs_l)){
  #hex_points <- spsample(borders, type = "hexagonal", cellsize = samps[i,1])
  #hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = cs_l[i])
  #fia_grid<-as.data.frame(poly.counts(aa.fia.pts, hex_grid))
  #cs_grid<-as.data.frame(poly.counts(aa.cs.pts, hex_grid))
  #res<-cs_l[i]
  #listid<-samps[i,2]
  #grid_pts<-cbind(fia_grid,cs_grid, rownames(fia_grid),res,listid)
  #names(grid_pts)<-c("fia","cs","id","res","listid")
  #grid_pts$diff<-round((grid_pts$fia/grid_pts$cs),2)
  #test<-hex_grid[c(grid_pts$id)]
  #my_list[[(length(my_list) + 1)]]<-as.data.frame(grid_pts)
  #my_list2[[length(my_list2) + 1]]<-hex_grid
#}
#library(plyr)
#my_df<-(ldply(my_list)[,])
#my_df<-my_df[my_df$diff>=0.9 & my_df$diff<=1.1,]
#my_df$sefia<-my_df$fia/my_df$res
#my_df$secs<-my_df$cs/my_df$res
#my_df<-unique(my_df)
#hist(my_df$sefia)
#my_df$hexidn<-as.numeric(sub("ID","",my_df$id))

#my_list3<-list()#extract polygons where sampling effort is even (fia vs. cs) and high (points per area) 
#for (i in 1:length(my_df$listid)){
  #p<-my_list2[[my_df$listid[i]]][my_df$hexidn[i]]
  #my_list3[[length(my_list3) + 1]]<-p
#}

#for (i in 1:length(my_list3)){
  #my_list3[[i]]@polygons[[1]]@ID<-paste0(my_df[i,3],"_",my_df[i,2])
#}

#evenhex<-SpatialPolygons(lapply(my_list3, function(x){x@polygons[[1]]}))
#evenhex<-gUnion(evenhex,evenhex)
#evenhex<-crop(evenhex,borders)
#plot(borders)
#plot(evenhex,add=TRUE)

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

#Plot predictors
#tiff('/~Desktop/myExpl_250m.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
#par(oma=c(0,0,1,2))
#plot(myExpl_250m[[1]],axes=F)
#plot(evenhex, lwd=2)
#dev.off()

#myExpl_df<-as.data.frame(myExpl_250m)# Check for multicollinearity
#M <- cor(na.omit(myExpl_df))
#tiff('~/Desktop/predcorrplot_evenhex_rem.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
#corrplot.mixed(M,tl.cex=0.5)#removed biovar6 and roads bc of multicollinearity issues
#dev.off()

#3. Rasterize response data to create presence and pseudoabsences----
#3A. FIA data
aa.cs.ras <- rasterize(x=aa.cs.pts, y=myExpl_250m[[1]], fun='count', background=0);
aa.cs.ras <- (aa.cs.ras*((myExpl_250m[[1]]*0)+1))>0
a2.cs.pts <- rasterToPoints(aa.cs.ras)
myRespName <- 'A_altissima_cs'
myResp <- a2.cs.pts[, 3] # the presence/absences data for our species
myResp[myResp==0] <- NA # setting 'true absences' to undefined
myRespXY <- a2.cs.pts[, c(1,2)] # the XY coordinates of species data
sum(na.omit(myResp))
myBiomodData.cs250m <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = myExpl_250m,
                                         resp.xy = myRespXY,
                                         resp.name = myRespName,
                                         PA.nb.rep = 1,
                                         PA.strategy = 'random',
                                         PA.nb.absences = sum(myResp, na.rm=T))
plot(myBiomodData.cs250m)
saveRDS(myBiomodData.cs250m, file = "myBiomodData.cs.rds")#3617 presences

#3Ai. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()
#3Aii. Computing the models
getwd()#outputs placed here
myBiomodModelOut.cs250m <- biomod2::BIOMOD_Modeling(myBiomodData.cs250m,
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
                                                 modeling.id=paste(myRespName,"FirstModeling.cs250m",sep=""))
saveRDS(myBiomodModelOut.cs250m, file = "myBiomodModelOut.cs250m.rds")
myBiomodModelEval.cs250m <- get_evaluations(myBiomodModelOut.cs250m) # get all models evaluation
dimnames(myBiomodModelEval.cs250m) # print the dimnames of this object
allmodels.cs250m.eval<-myBiomodModelEval.cs250m[c('TSS'),"Testing.data",,,] # print the eval scores of all selected models
allmodels.cs250m.eval<-as.data.frame(allmodels.cs250m.eval)
names(allmodels.cs250m.eval)<-"TSS"
models<-rownames(allmodels.cs250m.eval)
res<-250
db<-"cs"
allmodels<-cbind(models,as.data.frame(allmodels.cs250m.eval),res,db)
write.table(allmodels,"allmodels.csv", sep=",", row.names = FALSE,col.names = TRUE,append=TRUE)

vars_importance.cs250m <- as.data.frame(get_variables_importance(myBiomodModelOut.cs250m)) # print variable 
var_means<-rowMeans((vars_importance.cs250m),na.rm=TRUE)
res<-250
db<-"cs"
pred<-names(var_means)
var_means<-cbind(pred,as.data.frame(var_means),res,db)

#tiff('~/Desktop/tohexports/cs250mmeanvarimp.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
#barplot(colMeans(vars_imp_means_cs250m), ylab="Mean Variance Importance Score", xlab="Predictors")
#dev.off()
var_rank<-rank((vars_importance.cs250m)*-1,ties.method="average",na.last=NA)
var_rank<-cbind(pred,as.data.frame(var_rank),res,db)

write.table(var_means,"var_means.csv", row.names = FALSE, col.names = TRUE, append=TRUE, sep=",")#remove roads,canopy, and population
write.table(var_rank,"var_ranks.csv", row.names = FALSE, col.names = TRUE, append=TRUE, sep=",")#remove roads,canopy, and population

#3AiiiEnsembling the models
myBiomodEM.cs250m <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut.cs250m,
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

saveRDS(myBiomodEM.cs250m, file = "myBiomodEM.cs250m.rds")
ensmodeleval<-get_evaluations(myBiomodEM.cs250m) # get evaluation scores
ensmodeleval<-ensmodeleval[[2]]
stats<-row.names(ensmodeleval)
ensmodeleval<-cbind(stats,as.data.frame(ensmodeleval))
write.table(append = T, sep=",", ensmodeleval, "enseval.csv", row.names = FALSE, col.names = TRUE)

###3Aiv. projection over the globe under current conditions
myBiomodProj.cs250m <- BIOMOD_Projection(modeling.output = myBiomodModelOut.cs250m,
                                  new.env = myExpl_250m,
                                  proj.name = 'current.cs',
                                  selected.models = 'all',
                                  binary.meth = 'TSS',
                                  compress = 'xz',
                                  clamping.mask = F,
                                  output.format = '.grd')
saveRDS(myBiomodProj.cs250m, file = "myBiomodProj.cs250m.rds")
myCurrentProj.cs <- get_predictions(myBiomodProj.cs250m) # if you want to make custom plots, you can also get the projected map
myBiomodEF.cs250m <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM.cs250m,
                                         projection.output = myBiomodProj.cs250m)
saveRDS(myBiomodEF.cs250m, file = "myBiomodEF.cs250m.rds")
plot(myBiomodEF.cs250m)
pred.out.cs250m <- myBiomodEF.cs250m@proj@val[[1]]
pred.out.cs250m[]<-pred.out.cs250m[]/1000
plot(pred.out.cs250m)
writeRaster(pred.out.cs250m, filename = 'toh_ens_cs250m_wgs.tif', format="GTiff")
crs<-"+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" 
pred.out.cs250m.lcc<-(projectRaster(pred.out.cs250m,crs=crs))
writeRaster(pred.out.cs250m.lcc, filename = 'toh_ens_cs250m_lcc.tif', format="GTiff")

#3Av. output plot
rpts <- rasterToPoints(pred.out.cs250m.lcc)
rdf <- as.data.frame(rpts)
ggsdm <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
  geom_path(data=borders, aes(x=long, y=lat, group=group), col='white', lwd=1.1, alpha=.3) + 
  scale_fill_continuous(type='viridis') + theme_void() + theme(legend.position='none')
png(paste('pred.out.cs.250m.lcc', 
          gsub(':', '', substr(Sys.time(), 12, 19)), '.png', sep=''), 
    height=1080, width=2160); plot(ggsdm); dev.off()

