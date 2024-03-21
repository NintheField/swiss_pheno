# Reconstruct phenology based on hisotrical temperature data
# 1763 to 1960

# packages
library(ncdf4)
library(doParallel)
library(phenor)

# set variables
pv <- "maesh13d"
classes <- "1234_longseries_filled_NA_AUR"
burnin <- 4000
iter <- 36000
date <- "2024-03-13"

#### directories
base_dir <- file.path(
  "data/calibrations",
  paste0("bt_run_",date,"_class",classes,"_",burnin,"burnin_",iter,"iter")
)

recon_dir <- file.path(base_dir,pv)

### read rds file
files <- list.files(file.path(base_dir,pv), pattern = "model_comparison", full.names = T)
mod_nams <- paste0(pv,"_",stringr::str_sub(files,-12,-5))
mod_outs <- lapply(files, readRDS)
names(mod_outs) <- mod_nams

# extract temperature data for Liestal
liestal_coord <- c(2622208, 1260342)
calc <- F

if (calc) {

  tempfile <- list.files("../swiss_recon/temp_ens/", pattern = "_temp_EnKF_ens", full.names = T)
  ft <- nc_open(tempfile[1])
  N <- ncvar_get(ft, "N");E <- ncvar_get(ft, "E")
  nc_close(ft)

  dates <- seq(as.Date("1763-01-01"),as.Date("2020-12-31"), by = "day")
  lind <- c(which.min(abs(E - liestal_coord[1])),which.min(abs(N - liestal_coord[2])))

  temp_list <- array(NA, dim=c(length(tempfile),10,365))
  rownames(temp_list) <- stringr::str_sub(tempfile,48,51)

  # get closest grid cell of liestal
  for(ff in 1:(length(tempfile)-1)){

    ft1 <- nc_open(tempfile[ff])
    time1 <- as.Date(ncvar_get(ft1, "time"), origin = switch(grepl("EnKF",tempfile[ff]) + 1,"1900-01-01","1970-01-01"))
    varid <- switch(grepl("EnKF",tempfile[ff]) + 1,"TabsD","temp")
    ty1 <- ncvar_get(ft1,varid = varid,start = c(lind[1],lind[2],1,1), count =c(1,1,-1,-1))[,(length(time1) - 101):length(time1)]
    nc_close(ft1)

    ft2 <- nc_open(tempfile[ff + 1])
    varid <- switch(grepl("EnKF",tempfile[ff + 1]) + 1,"TabsD","temp")
    ty2 <- ncvar_get(ft2,varid = varid, start = c(lind[1],lind[2],1,1), count =c(1,1,-1,-1))[,1:263]
    nc_close(ft2)

    tynew <- cbind(ty1,ty2)
    rm(ty1,ty2)

    temp_list[ff + 1,,] <- tynew

  }

  saveRDS(temp_list, file = paste0("data-raw/temp/liestal_tempextraction_allmembers.rds"))

} else {
  temp_list <- readRDS(paste0("data-raw/temp/liestal_tempextraction_allmembers.rds"))
}

# get lat values for estimating the daylength
Ep <- (liestal_coord[1]-2600000)/1000000
Np <- (liestal_coord[2]-1200000)/1000000

latdegb_lt <- (16.9023892+3.238272*Np-0.270978*(Ep^2)-0.002528*(Np^2)-0.0447*(Ep^2)*Np-0.0140*(Np^3))*100/36
latvec <- as.vector(latdegb_lt)
Lidoy <- daylength(c(264:365,1:263),latvec)
doy <- c(-102:-1,1:263)

# dimensions of temp list
dd <- dim(temp_list)

registerDoParallel(cores = floor(detectCores()*0.2))
doymat <- array(NA, dim = c(dd[1],dd[2],length(mod_outs),length(mod_outs[[1]]$modelled)))

for(per in 1:(length(mod_outs))){

  # create list, ordered by year & ensemble
  temp_list_help <- split(temp_list,seq(dd[1]*dd[2]))

  # create list for input to phenor
  temp_help <- lapply(temp_list_help, function(temp) y <- list(Ti = matrix(temp, ncol = 1),
                                                               Li = matrix(Lidoy, ncol = 1),
                                                               doy = doy))

  mods <- names(mod_outs[[1]]$modelled)

  doyout <- array(sapply(1:length(mods), function(mm){

    # get model and parameters
    opt_params <- mod_outs[[per]]$modelled[[mm]]$parameters
    func <- get(mods[mm])

    doyvec <- unlist(foreach(a = 1:length(temp_help)) %dopar% {
      func(par = opt_params, data = temp_help[[a]])
    })

    return(doyvec)

  }), dim = c(dd[1],dd[2],length(mods)))

  doymat[,,per,] <- doyout
  doymat[doymat == 9999] <- NA

}

dimnames(doymat) <- list(years = rownames(temp_list),ens  = 1:dd[2], per = 1:length(mod_outs), mod = mods)

saveRDS(doymat, file = paste0(recon_dir,"/liestal_estimate_",pv,"_periods_model_ensemble_",classes,".rds"))
print("done")

