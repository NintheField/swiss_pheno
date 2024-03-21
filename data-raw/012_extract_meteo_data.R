# Extract temperature data from four closest grid cell
# Calculate mean temperature based on weights derived from topographic differences between series and grid cells

# packages
library(ncdf4)
#library(stats)
#library(geosphere)

# directories
indir_nc <- "../swiss_recon/temp/"

# functions
source("R/helpfuns.R")

# set variables
classes <- "1234_longseries_filled"
pv <- "mpica13d"
yrs <- 1951:2020
nyr <- length(yrs)

# read gridded data
ft <- nc_open(paste(indir_nc,"CH_temp_EnKF_1951-01-01-1951-12-31.nc",sep="")) # nonly for coordinates
lon <- ft$dim[[2]]$vals
lat <- ft$dim[[3]]$vals
nc_close(ft)

lonmat <- matrix(rep(lon,length(lat)),nrow=length(lon),ncol=length(lat))
latmat <- matrix(rep(lat,each = length(lon)),nrow=length(lon),ncol=length(lat))

# topo
file <- "../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc"
nc <- nc_open("../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc")
topo <- ncvar_get(nc,"height")

Ntopo <- ncvar_get(nc,"chy") + 1000000
Etopo <- ncvar_get(nc,"chx") + 2000000

# read pheno file
pheno_file <- paste0("data-raw/pheno_prep/pheno_dat_",pv,"_class",cllps(classes),".rds")
pheno_sel <- readRDS(pheno_file)

# extract station temperatures
meteo_dat <- array(NA,dim=c(nrow(pheno_sel),365))
selyrs <- unique(pheno_sel$reference_year)
tind_start <- 1:263

for (yr in selyrs){
  print(yr)

  # calculate the topo weights for every staation of reference year
  ahelp <- which(pheno_sel$reference_year== yr)
  pheno_sel_yr <- pheno_sel[ahelp,]
  nstat <- nrow(pheno_sel_yr)
  obsE <- pheno_sel_yr$obsE
  obsN <- pheno_sel_yr$obsN
  obsnam <- as.character(pheno_sel_yr$nat_abbr)
  obsalt <- pheno_sel_yr$alt

  sellon <- matrix(NA,nrow=nstat,ncol=4)
  sellat <- matrix(NA,nrow=nstat,ncol=4)
  topo_wd <- matrix(NA, nrow = nstat, ncol = 4)

  for (a in 1:nstat){
    dist <- (((obsE[a]-lonmat)^2)+(((obsN[a]-latmat))^2))^0.5
    sellon[a,] <- arrayInd(which(dist %in% dist[order(dist)[1:4]]),dim(dist))[1:4,1]
    sellat[a,] <- arrayInd(which(dist %in% dist[order(dist)[1:4]]),dim(dist))[1:4,2]
    topodiff <- abs(sapply(1:4, function(x) topo[sellon[a,x],sellat[a,x]]) - obsalt[a])
    topo_wd[a,] <- topodiff/sum(topodiff)
  }

  # get the data for the current year
  ft <- nc_open(paste(indir_nc,"CH_temp_TabsD_",yr,".nc",sep=""))
  ty <- ncvar_get(ft,varid="TabsD")
  nc_close(ft)

  tind_yrend <- (dim(ty)[3] - 101):dim(ty)[3]

  for (a in 1:nstat){
    temp_wd <- rowSums(sapply(1:4, function(x) ty[sellon[a,x],sellat[a,x],tind_start] * topo_wd[a,x] ))
    meteo_dat[ahelp[a],(length(tind_yrend)+1):365] <- round(temp_wd,2) # write the data from the current spring to doy 102 from current year
  }

  # get the data for the winter of the year before
  if(yr == 1961){
    ft <- nc_open(paste(indir_nc,"CH_temp_EnKF_",yr-1,"-01-01-",yr-1,"-12-31.nc",sep=""))
    ty <- ncvar_get(ft,varid="temp")
  } else{
    ft <- nc_open(paste(indir_nc,"CH_temp_TabsD_",yr-1,".nc",sep=""))
    ty <- ncvar_get(ft,varid="TabsD")
  }
  nc_close(ft)

  tind_yrend <- (dim(ty)[3] - 101):dim(ty)[3]

  for (a in 1:nstat){
    temp_wd <- rowSums(sapply(1:4, function(x) ty[sellon[a,x],sellat[a,x],tind_yrend] * topo_wd[a,x] ))
    meteo_dat[ahelp[a],1:length(tind_yrend)] <- round(temp_wd,2) ## write the data from the current winter to doy 1 of next year
  }
}

rownames(meteo_dat) <- as.character(pheno_sel$nat_abbr)

saveRDS(
  meteo_dat,
  file = paste0("data-raw/pheno_prep/meteo_dat_",pv,"_class",classes,".rds"),
  compress = "xz"
)

