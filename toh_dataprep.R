#set wd
setwd("~/Google Drive/My Drive/PhD/pops/slf/toh_exports/")
#load libraries
library(rinat)
library(rgbif)
library(rFIA)
library(sf)
library(data.table)
library(rgdal)
library(raster)
library(parallel)
library(dplyr)
cores<-detectCores()
options(timeout=1800)

#Setup reference area & compare to L48; setup species of interest
ref_area <- readOGR("/Volumes/GoogleDrive/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/larger_region_no_west_from_pseodo_mercator_in_pseudo_mercator.shp")
ref_area_t <- spTransform(ref_area, CRS("+proj=longlat +datum=WGS84 +no_defs"))
ref_area_xy<-as.vector(extent(ref_area_t))
l48<-readOGR("/Volumes/GoogleDrive/Shared drives/APHIS  Projects/shared resources/data/usa_boundaries/us_lower_48_states.shp")
l48_t<- spTransform(l48, CRS("+proj=longlat +datum=WGS84 +no_defs"))
plot(l48_t)
plot(ref_area_t, add=T)
sn <- "Ailanthus altissima"
#Import inat aa data, clean, convert to points
start_time <- Sys.time()
inat_out <- get_inat_obs(taxon_name = sn, maxresults = 10000, 
                         geo = TRUE, bounds = ref_area_t, quality = "research")
end_time <- Sys.time()
inat_b<- end_time - start_time

#Import rgbif data, clean, convert to points
start_time <- Sys.time()
gbif_out <- occ_data(scientificName = sn, hasCoordinate=TRUE,
                   geometry=c(ref_area_xy[1],ref_area_xy[3],ref_area_xy[2],ref_area_xy[4]), 
                   limit=100000, hasGeospatialIssue = FALSE)
end_time <- Sys.time()
gbif_b<- end_time - start_time
occ_aa<-as.data.table(gbif_out$data)
#Import EddMaps (downloaded 3/19), clean, convert to points
eddmaps_aa <- readOGR("/Volumes/GoogleDrive/My Drive/PhD/pops/slf/datasources/eddmaps/30797/observations.gpkg")
#Import rFIA, download separately, merge, clean, convert to points
start_time <- Sys.time()
fia <- getFIA(states = c('PA','NY','NJ','OH','WV','KY','VA','TN','NC','MD',
                            'DE'), nCores=(cores-1))
end_time <- Sys.time()
fia_b<- end_time - start_time
#Delete overlap between GBIF and inaturalist
gbif_out_data<-as.data.table(gbif_out[["data"]])
gbif_out_data<-gbif_out_data[institutionCode != "iNaturalist",]
gbif_out_data<-gbif_out_data[ ,`:=`(networkKeys = NULL)]
#Append unique ID to each TOH table and export .csv
gbif_out_data$index <- 1:nrow(gbif_out_data)
#write.csv2(gbif_out_data, "gbif_aa_all.csv")

inat_out$index<-1:nrow(inat_out)
#write.csv2(inat_out, "inat_aa_all.csv")

eddmaps_aa$index<-1:nrow(eddmaps_aa)
#write.csv2(eddmaps_aa, "eddmaps_aa_all.csv")

fia_tree<-as.data.table(fia["TREE"])
fia_plot<-as.data.table(fia["PLOT"])
fia_aa_tree<-fia_aa_tree[TREE.SPCD == 341,]
fia_aa<-merge(fia_aa_tree, fia_plot, by.x="TREE.PLT_CN", by.y="PLOT.CN", all.x=T)
fia_aa$index<-1:nrow(fia_aa)
fia_aa<-as.data.table(fia_aa)
#write.csv2(fia_aa, "fia_aa_all.csv")

#"clean" TOH tables to retain source database, record type, latitude, & longitude
fia_aa_c<-fia_aa[,.(uid=index,yr=TREE.INVYR,treeid=TREE.TREE,deadyr=TREE.MORTYR,
          lat=PLOT.LAT,lon=PLOT.LON, design=PLOT.DESIGNCD, db="fia")]
dups<-fia_aa_c[,.(yr=max(yr), uid=max(uid)), by=.(design,treeid)]
fia_aa_c<-merge(dups, fia_aa_c, by=c("treeid","design","yr","uid"))

gbif_aa_c<-gbif_out_data[basisOfRecord == "HUMAN_OBSERVATION" & year >1979,
                         .(uid=index,lat=decimalLatitude,
                           lon=decimalLongitude,yr=na.omit(year), db="gbif")]


inat_out<-as.data.table(inat_out)
inat_aa_c<-inat_out[coordinates_obscured == "false",
                    .(uid=index,lat=latitude,lon=longitude,yr=as.Date(datetime), db="inat")]
inat_aa_c$yr<-as.numeric(format(inat_aa_c$yr,format="%Y"))

eddmaps_aa<-as.data.table(eddmaps_aa)
eddmaps_aa_c<-eddmaps_aa[OccStatus =="Detected" & Voucher !=1,
                                .(uid=index,yr=as.Date(ObsDate, format="%m-%d-%y"),lat=coords.x2,
                                  lon=coords.x1,db="eddmaps")]
eddmaps_aa_c$yr<-as.numeric(format(eddmaps_aa_c$yr,format="%Y"))
eddmaps_aa_c[yr>1980,]

#Create & export: FIA, citizen-science, & both data tables; set aside 20% for validation
#FIA
testing<-as.numeric(sample(row.names(fia_aa_c),size=round((nrow(fia_aa_c)*0.8),0)))
fia_ib<-fia_aa_c[testing,.(uid,db,lat,lon)]
fia_oob<-fia_aa_c[!testing,.(uid,db,lat,lon)]

#write.csv2(fia_ib,"fia_ib.csv")
#write.csv2(fia_oob,"fia_oob.csv")
#Citizen-science (GBIF, inaturalist, & EddMaps)

#GBIF
testing<-as.numeric(sample(row.names(gbif_aa_c),size=round((nrow(gbif_aa_c)*0.8),0)))
gbif_ib<-gbif_aa_c[testing,.(uid,db,lat,lon)]
gbif_oob<-gbif_aa_c[!testing,.(uid,db,lat,lon)]
#write.csv2(gbif_ib,"gbif_ib.csv")
#write.csv2(gbif_oob,"gbif_oob.csv")
#iNaturalist
testing<-as.numeric(sample(row.names(inat_aa_c),size=round((nrow(inat_aa_c)*0.8),0)))
inat_ib<-inat_aa_c[testing,.(uid,db,lat,lon)]
inat_oob<-inat_aa_c[!testing,.(uid,db,lat,lon)]
#write.csv2(inat_ib,"inat_ib.csv")
#write.csv2(inat_oob,"inat_oob.csv")
#EddMaps
testing<-as.numeric(sample(row.names(eddmaps_aa_c),size=round((nrow(eddmaps_aa_c)*0.8),0)))
eddmaps_ib<-eddmaps_aa_c[testing,.(uid,db,lat,lon)]
eddmaps_oob<-eddmaps_aa_c[!testing,.(uid,db,lat,lon)]
#write.csv2(eddmaps_ib,"eddmaps_ib.csv")
#write.csv2(eddmaps_oob,"eddmaps_oob.csv")
#Citizen-science data merged and exported
cs_ib<-rbind(eddmaps_ib,inat_ib,gbif_ib)
cs_oob<-rbind(eddmaps_oob,inat_oob,gbif_oob)
#write.csv2(cs_ib,"cs_ib.csv")
#write.csv2(cs_oob,"cs_oob.csv")
#Both (merge FIA & citizen-science)
all_ib<-rbind(cs_ib,fia_ib)
all_oob<-rbind(cs_oob,fia_oob)
#write.csv2(cs_ib,"cs_ib.csv")
#write.csv2(cs_oob,"cs_oob.csv")
