find_edges <- function(map,n=3)   {

  mapold <- map
  neighbours <- matrix(0,nrow=nrow(map),ncol=ncol(map))

  mat <- matrix(NA,nrow=n,ncol=ncol(map))
  map <- rbind(mat,map,mat)

  mat <- matrix(NA,nrow=nrow(map),ncol=n)
  map <- cbind(mat,map,mat)

  N <- nrow(map); M <- ncol(map)

  arr <- array(dim=c(N-2*n,M-2*n,4*n))
  for (i in 1:n)     arr[,,(4*(i-1)+1):(4*i)] <- array( c( map[(n+1-i):(N-n-i),(n+1):(M-n)], map[(n+1+i):(N-n+i),(n+1):(M-n)],
                                                           map[(n+1):(N-n),(n+1-i):(M-n-i)], map[(n+1):(N-n),(n+1+i):(M-n+i)]), dim=c(N-2*n,M-2*n,4) )

  neighbours <- apply(arr,c(1,2),mean,na.rm=TRUE)

        edges = neighbours < 1 & mapold;
        return(edges)
}


nlev=1:length(lev)

# 47 deciduous forest edge
ndf <- nlev[lev=='Deciduous forest']
uds_mat[find_edges(uds_mat==ndf)==1] = 47;
lev[47]<-"Deciduous forest edge"


# 48 olive groves edge
ndf <- nlev[lev=='olive groves']
uds_mat[find_edges(uds_mat==ndf)==1] = 48;
lev[48]="olive groves edge"


# 49 Orchards and small fruit farms edge
ndf <- nlev[lev=='Orchards and small fruit farms']
uds_mat[find_edges(uds_mat==ndf)==1] = 49;
lev[49]="Orchards and small fruit farms edge"


# 50 Early growth forests (reforestation) edge
ndf <- nlev[lev=='Early growth forests (reforestation)']
uds_mat[find_edges(uds_mat==ndf)==1] = 50;
lev[50]="Early growth forests (reforestation) edge"


# 51 Mixed coniferous and decidous forest edge
ndf <- nlev[lev=='Mixed coniferous and decidous forest']
uds_mat[find_edges(uds_mat==ndf)==1] = 51;
lev[51]="Mixed coniferous and decidous forest edge"


# 52 coniferous forest edge
ndf <- nlev[lev=='coniferous forest']
uds_mat[find_edges(uds_mat==ndf)==1] = 52;
lev[52]="coniferous forest edge"


# 53 Agricultural production units edge
ndf <- nlev[lev=='Agricultural production units']
uds_mat[find_edges(uds_mat==ndf)==1] = 53;
lev[53]="Agricultural production units edge"


# 54 Other perennial crops edge
ndf <- nlev[lev=='Other perennial crops']
uds_mat[find_edges(uds_mat==ndf)==1] = 54;
lev[54]="Other perennial crops edge"


# 55 vineyards edge
ndf <- nlev[lev=='vineyards']
uds_mat[find_edges(uds_mat==ndf)==1] = 55;
lev[55]="vineyards edge"


# 56 Cultivation and complex systems edge
ndf <- nlev[lev=='Cultivation and complex systems']
uds_mat[find_edges(uds_mat==ndf)==1] = 56;
lev[56]="Cultivation and complex systems edge"



