rm(list=ls()); graphics.off()

# Load packages
library(rgeos)
library(rgdal)
library(raster)
library(rts)
library(maptools)   # for shp files
library(Matrix)     # for sparse matrices
library(pracma)     # for repmat
library(ggplot2)
library(plyr)

dx <- 20

out_name <- paste("output",dx,sep="")

suit_name <- paste("data/HS_CVdx",dx,".tif",sep="")
presen_path <- paste("data/pres",dx,".tif",sep="")
domain_path <- "data/boundaryPA.shp"
param_path  <- "data/parameters.csv"
uds_path    <- "data/land_cover.shp"

# Read inputs files
param    <- read.csv(param_path,sep=';',header=FALSE)
pres     <- raster(presen_path)
suppressWarnings( domain   <- readShapeSpatial(domain_path,proj4string=CRS("+proj=longlat")) )

# Parameters of the model
par_values <- as.character(param$V2)
D      <- as.numeric(par_values[1])  
r      <- as.numeric(par_values[2])  
c<-1
T <- 18

date0 <- as.Date("2012/7/1")

k <-1; q<- 2;
m <- 2*q-1; #alpha <- c*m/(B^m)

# Numerical grid

dxkm <- dx*1e-3

yr<-2012

rasterOptions(chunksize=1e+4, maxmemory=1e+4, progress='text', overwrite=TRUE)
#pres_num <- raster::aggregate(pres, fact=dx/xres(pres), fun=max, extend=FALSE)  #!!!!! max no mean
#pres_num <- raster('data/pres10.tif')

# Suitability

   HSI_rast <- raster(suit_name)
   HSI <- getValues(HSI_rast,format="matrix")
   HSI <- t(HSI)
   HSI <- HSI[,ncol(HSI):1]
   HSI <- c(HSI)

pres_num <- pres * HSI_rast
#plot(pres_num)
#lines(domain,col=2)
#title(paste(yr))


# Tolerance settings
tolsolve <- 1e-15


# Time settings
dt <- 0.05
int <- 1/dt



# Calculate some values
Nx <- dim(pres_num)[2]
Ny <- dim(pres_num)[1]
n <- Nx*Ny
Nt <- round(T/dt) +1
t <- ((0:(Nt-1)))*dt
#elo <- exp(-delta*dt)


### Assemby matrices
load(paste("L",dx,".RData",sep=""))

# Matrix L (without 1/dxkm^2 factor)
L = bandSparse(n,n, k=c(0,Nx),diagonals=list(4+numeric(n),-1+numeric(Nx*(Ny-1))), symmetric = T)
# S
L[1,1]=3; L[1,2]=-3/2; L[Nx,Nx]=6; L[Nx,Nx-1]=-3
L[cbind(2:(Nx-1),1:(Nx-2))]=-1;
L[cbind(2:(Nx-1),3:Nx)]=-1;
#Z
L[(Ny-1)*Nx+1,(Ny-1)*Nx+1]=6; L[n,n]=3
L[(Ny-1)*Nx+1,(Ny-1)*Nx+2]=-3; L[n,n-1]=-3/2
L[cbind(((Ny-1)*Nx+2):(Nx*Ny-1),((Ny-1)*Nx+1):(n-2))]=-1
L[cbind(((Ny-1)*Nx+2):(Nx*Ny-1),((Ny-1)*Nx+3):(n))]=-1
#T
L[1,Nx+1]=-3/2; L[Nx,2*Nx]=-3
L[cbind(2:(Nx-1),(Nx+2):(2*Nx-1))]=-2
#Y
L[(Ny-1)*Nx+1,(Ny-2)*Nx+1]=-3; L[n,(Ny-1)*Nx]=-3/2
L[cbind( ((Ny-1)*Nx+2):(n-1), ((Ny-2)*Nx+2):((Ny-1)*Nx-1))]=-2
#X
ind1 <- ind2 <- matrix(nrow=0,ncol=2)
for (j in 1:(Ny-2)) {
    ind1 <- rbind(ind1, cbind((j*Nx+2):((j+1)*Nx-1),(j*Nx+1):((j+1)*Nx-2)),cbind((j*Nx+2):((j+1)*Nx-1),(j*Nx+3):((j+1)*Nx)) )
    ind2 <- rbind(ind2, c(j*Nx+1,j*Nx+2), c((j+1)*Nx,(j+1)*Nx-1))
}
L[ind1]=-1; L[ind2]=-2

# save(L, file=paste("L",dx,".RData",sep=""))


# Matrix B
Bmat <- bandSparse(n,n, k=0,diagonals=list(1+numeric(n))) + dt/(dxkm^2)*D*L

# Initial conditions
U0 <- getValues(pres_num,format="matrix")
U0 <- t(U0)
U0 <- U0[,ncol(U0):1]
u0 <- c(U0)

uout <- repmat(matrix(u0,nrow=n),1,T+1)


    # Forward
    u = u0

    for (nt in 2:Nt) {

    print(nt)

        Fun =  r*HSI*u-r*u^2/k #- mu*u*Evec/(1+h*mu*u)
        u = u + dt*Fun

        u <- solve(Bmat,u,sparse=TRUE,tol=tolsolve)

        if ( mod(nt-1,int) ==0 )  uout[,(nt-1)/int+1] <- as.numeric(u)
         
    }


### Output

dates <- seq(date0, by = "year", length.out = T+1)

dir.create(paste("data/",out_name,sep="") )

K<-5e+5
A<-dx^2*1e-6
KA<-K*A


for (s in 1:ncol(uout)) {
    us <- uout[,s]
    Us <- matrix(us,ncol=Ny)
    Us <- Us[,ncol(Us):1]
    Us <- t(Us)
    dens <- setValues(pres_num, Us, format='matrix')
    pres <- dens*KA
    disp(yr)
    writeRaster(pres,paste("data/",out_name,"/","presence_dx",dx,"_",yr,sep=""),format="GTiff", overwrite=TRUE)
    writeRaster(dens,paste("data/",out_name,"/","density_dx",dx,"_",yr,sep=""),format="GTiff", overwrite=TRUE)
    writeRaster(dens/HSI_rast,paste("data/",out_name,"/","densityHSIratio_dx",dx,"_",yr,sep=""),format="GTiff", overwrite=TRUE)
    windows(); plot(dens); title(yr)
    yr<-yr+1
    #ab_rast <- stack(ab_rast, dens)
}

