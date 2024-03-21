# Select longest series depending on pheno phase,class, below certain altitude
# Calculate climatic bins for all series of certain class
# Fill missing values with observations from the same climatic bin

#### Packages
library(ncdf4)
library(lubridate)

# functions
source("R/helpfuns.R")

# directories
indir_nc <- "../swiss_recon/temp/"

# set variables
yrs <- 1991:2020
nyr <- length(yrs)
classes <- c(1,2,3,4)
topo_sel <- 1000
pv <- "mpica13d"
periods <- list(c(1961,1980),c(1971,1990),c(1981,2000),c(1991,2010),c(2001,2020))

# read station list
station_list <- read.csv("data-raw/pheno_net/data/ch.meteoschweiz.messnetz-phaenologie_de.csv",header = T, sep = ",", allowEscapes = T)
row.names(station_list) <- as.character(station_list$Abk.)

# read pheno data
pheno <- read.table("data-raw/pheno_net/phaeno_previous.csv",header=T, sep =";")

# read quality classification
classification <- read.table("data-raw/pheno_net/dataD4.csv",header=T, sep =",")

# read qc flags
qc <- read.csv("data-raw/pheno_net/PhenoClass/QCflags_minimum_Final_Juni2017.csv",header=T, sep =",")

# phenovals and param_id
phenovals <- c("mfags13d", "mprua65d", "maesh13d","mlard13d","mpica13d")
param_id <- c(606,662,601,634,637)
names(param_id) <- phenovals

# gridded data
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

# select series based on param_id and quality classes
vval_stns <- as.character(classification$station[which(classification$OLD_SHORT_TX == pv & classification$Rank %in% classes)])
if(!is.na(topo_sel)){vval_stns  <- station_list$Abk.[station_list$Stationshoehe < topo_sel & station_list$Abk.%in% vval_stns]}

pheno_sel <- pheno[pheno$param_id == param_id[pv], ]
vval_stns <- vval_stns[vval_stns %in% unique(pheno_sel$nat_abbr)]
pheno_sel <- pheno_sel[as.character(pheno_sel$nat_abbr) %in% vval_stns,]

# add qc flags ## most of them are already remove though
qc_pv <- qc[qc$OLD_SHORT_TX == pv & qc$FinalFlagCountUnique == 1 & qc$station %in% vval_stns,c("station","year","doy","value","FinalFlagCountUnique")]

for(ii in 1:nrow(qc_pv)){
  ind <- which(pheno_sel$nat_abbr == qc_pv[ii,"station"] & pheno_sel$reference_year == qc_pv[ii,"year"])
  if(length(ind) > 0){
    print(pheno_sel[ind,])
    pheno_sel[ind,c("value","doy")] <- NA
  }
}

# select series based on start and end years
start_stns <- sapply(vval_stns, function(x){any(pheno_sel$reference_year[pheno_sel$nat_abbr == x] < 1966 & pheno_sel$reference_year[pheno_sel$nat_abbr == x] > 1960)})
end_stns <- sapply(vval_stns, function(x){any(pheno_sel$reference_year[pheno_sel$nat_abbr == x] > 2010)})

longseries <- vval_stns[(start_stns + end_stns) == 2]
pheno_sel_long <- pheno_sel[pheno_sel$nat_abbr %in% longseries,] ## from these only select long series
pheno_sel_long <- pheno_sel_long[pheno_sel_long$reference_year > 1960 & pheno_sel_long$reference_year < 2021,]

# create complete data frame
sites <- unique(pheno_sel_long$nat_abbr)
pheno_sel_full <- data.frame(matrix(NA, nrow= length(sites) * length(1961:2020), ncol = ncol(pheno_sel_long)))
colnames(pheno_sel_full) <- colnames(pheno_sel_long)
pheno_sel_full$reference_year <- rep(1961:2020,length(sites))
pheno_sel_full$nat_abbr <- rep(sites,each = length(1961:2020))
pheno_sel_full$param_id <- pheno_sel_long$param_id[1]

for(ss in 1:length(sites)){
  sind <- pheno_sel_long$nat_abbr == sites[ss]
  yind <- which(pheno_sel_full$reference_year %in% pheno_sel_long$reference_year[sind] & pheno_sel_full$nat_abbr == sites[ss])
  pheno_sel_full[yind,c("value","doy")] <- pheno_sel_long[sind,c("value","doy")]
}

# get mean climate of all vval series from selected classes, and put in bins based on their mean temperature
station_list_sel <- station_list[as.character(station_list$Abk.) %in% vval_stns,1:10]
rownames(station_list_sel) <- station_list_sel$Abk.

nstat <- nrow(station_list_sel)
obslon <- station_list_sel$Laengengrad
obslat <- station_list_sel$Breitengrad
obsE <- station_list_sel$KoordinatenE
obsN <- station_list_sel$KoordinatenN
obsnam <- as.character(station_list_sel$Abk.)
obsalt <- station_list_sel$Stationshoehe

sellon <- matrix(NA,nrow=nstat,ncol=4)
sellat <- matrix(NA,nrow=nstat,ncol=4)

topo_wd <- matrix(NA, nrow = nstat, ncol = 4)

# calculate topographic weight for series
for (a in 1:nstat){
  dist <- (((obsE[a]-lonmat)^2)+(((obsN[a]-latmat))^2))^0.5
  sellon[a,] <- arrayInd(which(dist %in% dist[order(dist)[1:4]]),dim(dist))[1:4,1]
  sellat[a,] <- arrayInd(which(dist %in% dist[order(dist)[1:4]]),dim(dist))[1:4,2]
  topodiff <- abs(sapply(1:4, function(x) topo[sellon[a,x],sellat[a,x]]) - obsalt[a])
  topo_wd[a,] <- topodiff/sum(topodiff)
}

# extract station temperatures
meteo_dat <- array(NA,dim=c(nstat,nyr,263))
startyr <- 1991

for (yr in startyr:2020){
  print(yr)

  ft <- nc_open(paste(indir_nc,"CH_temp_TabsD_",yr,".nc",sep=""))
  ty <- ncvar_get(ft,varid="TabsD")
  tind_start <- 1:263

  for (a in 1:nstat){
    temp_wd <- rowSums(sapply(1:4, function(x) ty[sellon[a,x],sellat[a,x],tind_start] * topo_wd[a,x] ))
    meteo_dat[a,yr - startyr + 1,] <- round(temp_wd,2) # write the data from the current spring to doy 102 from current year
  }
  nc_close(ft)

}
rownames(meteo_dat) <- obsnam

# calcualte same climatic bins
dates <- seq(as.Date("1991-01-01"),as.Date("1991-10-01"), by = "day")[1:263]
reldays <- month(dates) %in% c(2:4)
meteo_clim <- round(apply(meteo_dat[,,reldays],c(1),mean),2)
bins <- as.integer(cut(meteo_clim, breaks = seq(0,10,1)))
names(bins) <- obsnam

# make samples based these bins
pheno_sel$climbin <- bins[pheno_sel$nat_abbr]

# check for series with missing data
naind <- which(is.na(pheno_sel_full$doy))
pheno_sel_miss <- pheno_sel_full[naind,]

# if available fill missing observations with other series
newsites <- sites
for(nr in 1:nrow(pheno_sel_miss)){

  missbin <- bins[pheno_sel_miss[nr,"nat_abbr"]]
  binstns <- names(bins[bins == missbin])
  missyr <- pheno_sel_miss[nr,"reference_year"]

  fsel <- pheno_sel[pheno_sel$reference_year == missyr & pheno_sel$nat_abbr %in% binstns & !pheno_sel$nat_abbr %in% longseries,]
  if(!nrow(fsel) <1){
    repl_class <- classification[which(classification$station %in% fsel$nat_abbr & classification$OLD_SHORT_TX == pv),c("station","Rank")]
    rownames(repl_class) <- repl_class$station
    fsel$class <- repl_class[fsel$nat_abbr,"Rank"]

    # check whether a station already has been used
    i <- 1
    while(fsel[order(fsel$class),][i,"nat_abbr"] %in% newsites){
      if(i == nrow(fsel)){
        break
      } else{
        i <- i + 1
        print(i)
      }
    }

    repl_abbr <- fsel[order(fsel$class),][i,"nat_abbr"]

    pheno_sel_full[naind[nr],] <- pheno_sel[pheno_sel$reference_year == missyr & pheno_sel$nat_abbr == repl_abbr,1:5]
    newsites <- unique(pheno_sel_full$nat_abbr)
  }
}

# how many nas per period are still in there?
sapply(periods, function(x) sum(is.na(pheno_sel_full[pheno_sel_full$reference_year %in% x[1]:x[2],"doy"])))
sapply(periods, function(x) sum(!is.na(pheno_sel_full[pheno_sel_full$reference_year %in% x[1]:x[2],"doy"])))
sapply(periods, function(x) table(pheno_sel_full[pheno_sel_full$reference_year %in% x[1]:x[2],"nat_abbr"]))

# create a new station list with all variables
allseries <- unique(pheno_sel_full$nat_abbr)
station_list_full <- station_list[as.character(station_list$Abk.) %in% allseries,1:10]
station_list_full <- station_list_full[allseries,] # reorder station list according to vvals

# get the pheno data
stnidx <- match(pheno_sel_full$nat_abbr, station_list_full$Abk.)

lat <- station_list_full$Breitengrad[stnidx]
lon <- station_list_full$Laengengrad[stnidx]
N <- station_list_full$KoordinatenN[stnidx]
E <- station_list_full$KoordinatenE[stnidx]
alt <- station_list_full$Stationshoehe[stnidx]
pheno_dat <- data.frame(pheno_sel_full[,1:5], obsLon= lon, obsLat= lat, obsE= E, obsN=N, alt =alt)

# save as compressed RDS
saveRDS(
  pheno_dat,
  file = paste0("data-raw/pheno_prep/pheno_dat_",pv,"_class",cllps(classes),"_longseries_filled.rds"),
  compress = "xz"
)
