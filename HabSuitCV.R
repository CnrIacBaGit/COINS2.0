### Load packages ###

library(raster)
library(maptools)
rm(list=ls()); graphics.off()

### Read files ###
rasterOptions(chunksize=1e+4, maxmemory=1e+4, progress='text', overwrite=TRUE)
pres     <- raster("data/initial_density.tif")
suppressWarnings(uds <- readShapeSpatial("data/land_cover.shp",proj4string=CRS("+proj=longlat"))  )

dx <- 20
pres <- raster::aggregate(pres, fact=dx/xres(pres), fun=max, extend=FALSE)
plot(pres)
writeRaster(pres,paste("data/pres",dx,sep=""),format="GTiff", overwrite=TRUE)

pres_thres <- 0.1
lev <- levels(uds@data$DESC_ENG)
uds_rast <- rasterize(uds,pres,'DESC_ENG',fun="first")
writeRaster(uds_rast,paste("data/uds_rast",dx,sep=""),format="GTiff", overwrite=TRUE)

ext <- pres@extent
brk <- seq(ext@xmin,ext@xmax,length.out=6)
pres_list <- list(pres,pres,pres,pres,pres)
uds_list <- list(uds_rast,uds_rast,uds_rast,uds_rast,uds_rast)
for (i in 1:5) {
    ext@xmin <- brk[i]
    ext@xmax <- brk[i+1]
    pres_list[[i]] <- crop(pres,ext)
    uds_list[[i]] <- crop(uds_rast,ext)
}
rast <- list(pres,pres,pres,pres,pres)
rast[[1]] <- merge(pres_list[[2]],pres_list[[3]],pres_list[[4]],pres_list[[5]])
rast[[2]] <- merge(pres_list[[1]],pres_list[[3]],pres_list[[4]],pres_list[[5]])
rast[[3]] <- merge(pres_list[[1]],pres_list[[3]],pres_list[[4]],pres_list[[5]])
rast[[4]] <- merge(pres_list[[1]],pres_list[[2]],pres_list[[3]],pres_list[[5]])
rast[[5]] <- merge(pres_list[[1]],pres_list[[2]],pres_list[[3]],pres_list[[4]])

uds <- list(uds_rast,uds_rast,uds_rast,uds_rast,uds_rast)
uds[[1]] <- merge(uds_list[[2]],uds_list[[3]],uds_list[[4]],uds_list[[5]])
uds[[2]] <- merge(uds_list[[1]],uds_list[[3]],uds_list[[4]],uds_list[[5]])
uds[[3]] <- merge(uds_list[[1]],uds_list[[3]],uds_list[[4]],uds_list[[5]])
uds[[4]] <- merge(uds_list[[1]],uds_list[[2]],uds_list[[3]],uds_list[[5]])
uds[[5]] <- merge(uds_list[[1]],uds_list[[2]],uds_list[[3]],uds_list[[4]])


for (j in 1:5) {

  pres_par <- rast[[j]]
  uds_rast<-uds[[j]]

  pres_mat <- getValues(pres_par, format="matrix")
  pres_mat[pres_mat>pres_thres]  <- 1
  pres_mat[pres_mat<=pres_thres] <- 0
  uds_mat  <- getValues(uds_rast, format="matrix")
  Nrow=dim(pres_mat)[1]
  Ncol=dim(pres_mat)[2]

  source('edges_n.r')

  uds_rast_ext <- uds_rast
  uds_rast_ext[] <- uds_mat
  plot(uds_rast_ext)
  #writeRaster(uds_rast_ext,"data/LC_ext",format="GTiff", overwrite=TRUE)

  ind_pres <- which(pres_mat>0,arr.in=TRUE)
  ind_pres <- ind_pres[ind_pres[,1]>1 & ind_pres[,2]>1 & ind_pres[,1]<Nrow & ind_pres[,2]<Ncol,]
  #pres_mat[ind_pres]

  score <- numeric(length(lev))
  for (i in 1:nrow(ind_pres))  {
      ri <- ind_pres[i,1]
      ci <- ind_pres[i,2]
      score[uds_mat[ri,ci-1]] <- score[uds_mat[ri,ci-1]]+ ( pres_mat[ri,ci-1]==0 & !is.na(uds_mat[ri,ci-1]) )
      score[uds_mat[ri,ci+1]] <- score[uds_mat[ri,ci+1]]+ ( pres_mat[ri,ci+1]==0 & !is.na(uds_mat[ri,ci+1]) )
      score[uds_mat[ri-1,ci]] <- score[uds_mat[ri-1,ci]]+ ( pres_mat[ri-1,ci]==0 & !is.na(uds_mat[ri-1,ci]) )
      score[uds_mat[ri+1,ci]] <- score[uds_mat[ri+1,ci]]+ ( pres_mat[ri+1,ci]==0 & !is.na(uds_mat[ri+1,ci]) )
  }

  freq_norm <- score/max(score,na.rm=T)

  aree <- numeric(length(lev))
  for (l in 1:length(lev))  aree[l] <- sum(uds_mat==l,na.rm=TRUE)

  dens <- numeric(length(lev))
  dens[aree>0] <- score[aree>0]/aree[aree>0]
  dens_norm <- dens/max(dens,na.rm=T)

  tab <- data.frame(lev,freq=score,freq_norm,dens,dens_norm,area=aree)
  write.table(tab,paste("data/scoreCV",j,".csv",sep=""),row.name=FALSE,sep=";")

  rho_freq <- rho_dens <- uds_mat
  for (l in 1:length(lev))   { rho_freq[uds_mat==l] <- freq_norm[l]
                               rho_dens[uds_mat==l] <- dens_norm[l]
  }
  rho_freq[is.na(rho_freq)]<- 1
  rho_dens[is.na(rho_dens)]<- 1

  rhof <- rhod <- uds_rast
  rhof[] <-  rho_freq
  rhod[] <-  rho_dens

  writeRaster(rhof,paste("data/habitat_suitability_frequencyCV",j,"dx",dx,sep=""),format="GTiff", overwrite=TRUE)
  writeRaster(rhod,paste("data/habitat_suitability_densityCV",j,"dx",dx,sep=""),format="GTiff", overwrite=TRUE)
}


hs1 <- read.table(paste("data/scoreCV",1,".csv",sep=""),header=TRUE,sep=";") [,4]
hs2 <- read.table(paste("data/scoreCV",2,".csv",sep=""),header=TRUE,sep=";") [,4]
hs3 <- read.table(paste("data/scoreCV",3,".csv",sep=""),header=TRUE,sep=";") [,4]
hs4 <- read.table(paste("data/scoreCV",4,".csv",sep=""),header=TRUE,sep=";") [,4]
hs5 <- read.table(paste("data/scoreCV",5,".csv",sep=""),header=TRUE,sep=";") [,4]

lev <- read.table(paste("data/scoreCV",1,".csv",sep=""),header=TRUE,sep=";") [,1]
#aree <- read.table(paste("data/scoreCV",1,".csv",sep=""),header=TRUE,sep=";") [,6]

tab <- data.frame(hs1,hs2,hs3,hs4,hs5)
rownames(tab) <- lev
nlev <- length(lev)

HS <- rowMeans(tab, na.rm=T)
cbind(sort(HS))

sdHS <- numeric(nlev)
names(sdHS)<- lev
for (i in 1:nlev) sdHS[i] <- sd(tab[i,],na.rm=T)
cbind(sort(sdHS))

tabHS <- data.frame(mean=HS,sd=sdHS)
tabHS <- cbind(tab,tabHS)
write.table(tabHS,paste("data/HS_CV.csv",sep=""),sep=";",col.names = NA, qmethod = "double")


uds_rast <- raster(paste("data/uds_rast",dx,".tif",sep=""))
uds_mat  <- getValues(uds_rast, format="matrix")

rho_mat <- uds_mat
for (l in 1:nlev)   rho_mat[uds_mat==l] <- HS[l]
rho_mat[is.na(rho_mat)]<- 0

rho_rast <- uds_rast
rho_rast[] <-  rho_mat

writeRaster(rho_rast,paste("data/HS_CV"),format="GTiff", overwrite=TRUE)